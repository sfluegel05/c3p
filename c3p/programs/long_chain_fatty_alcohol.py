"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: Long-chain fatty alcohol
Definition: A fatty alcohol with a carbon chain length ranging from C13 to C22.
A valid long-chain fatty alcohol should contain an alkyl chain (of 13–22 carbons in a contiguous, acyclic path)
and at least one –OH substituent (satisfying a simple alcohol pattern) attached to one of the carbons in that chain.
Note: This is a heuristic method.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    We require that the molecule has a contiguous chain of carbon atoms (simple path) of length 13 to 22,
    with at least one carbon in that chain attached to an -OH group (alcohol).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule qualifies as a long-chain fatty alcohol, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string and add explicit hydrogens so we can see OH groups.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    molH = Chem.AddHs(mol)

    # Build a dictionary for all carbon atoms (atomic number 6) in molH.
    # Keys are atom indices; values will store whether that carbon is attached to an -OH group.
    carbon_atoms = {}
    for atom in molH.GetAtoms():
        if atom.GetAtomicNum() == 6:
            carbon_atoms[atom.GetIdx()] = False

    # Determine which carbon atoms carry an -OH substituent.
    # We consider an -OH group when an oxygen (atomic num 8) is attached, and the oxygen
    # has exactly one neighbor that is hydrogen (atomic num 1) and one neighbor that is carbon.
    for atom in molH.GetAtoms():
        if atom.GetAtomicNum() == 6:
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    # Check that oxygen is not, say, an ether: we expect O to have exactly 2 neighbors,
                    # one being hydrogen and one being the carbon we are looking at.
                    if len(nbr.GetNeighbors()) == 2:
                        hasH = any(n.GetAtomicNum() == 1 for n in nbr.GetNeighbors())
                        if hasH:
                            # Mark the carbon as having an attached -OH group.
                            carbon_atoms[atom.GetIdx()] = True
                            # (Once marked, one attached alcohol group is enough.)
                            break

    # Create a graph that includes only the carbon atoms.
    # The graph is represented as a dictionary: {atom_idx: [neighbor_indices, ...], ...}
    carbon_graph = {idx: [] for idx in carbon_atoms.keys()}
    for bond in molH.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            i1 = a1.GetIdx()
            i2 = a2.GetIdx()
            # Only add if both atoms are in our carbon_atoms dictionary.
            if i1 in carbon_graph and i2 in carbon_graph:
                carbon_graph[i1].append(i2)
                carbon_graph[i2].append(i1)

    # We now search for any simple (acyclic) path in the carbon_graph that has length in [13,22]
    # and for which at least one carbon in that path has an attached -OH group.
    found = False
    found_path_length = 0
    reason = ""
    max_path_length_found = 0  # track longest chain we see

    # Use DFS recursion to enumerate paths. Because the graphs are small,
    # a brute-force DFS is acceptable.
    def dfs(current, visited, path):
        nonlocal found, found_path_length, max_path_length_found
        path.append(current)
        visited.add(current)
        # Update longest chain length seen
        if len(path) > max_path_length_found:
            max_path_length_found = len(path)
        # Check if current path length is within desired range.
        if 13 <= len(path) <= 22:
            # Check if any carbon in this path has an -OH substituent.
            if any(carbon_atoms[atom_idx] for atom_idx in path):
                found = True
                found_path_length = len(path)
                # Once one chain qualifies, we can stop search.
                # We break out by clearing visited set (or simply return).
                # We do not call further recursion.
                visited.remove(current)
                path.pop()
                return
        # If path is already 22, do not extend further.
        if len(path) >= 22:
            visited.remove(current)
            path.pop()
            return

        # Continue DFS for unvisited neighbors.
        for neighbor in carbon_graph.get(current, []):
            if neighbor not in visited:
                dfs(neighbor, visited, path)
                if found:
                    # If a matching chain has been found, no need to search more.
                    return
        visited.remove(current)
        path.pop()

    # Start DFS from each carbon atom.
    for node in list(carbon_graph.keys()):
        if found:
            break
        dfs(node, set(), [])

    if found:
        return True, f"Found a carbon chain of length {found_path_length} with an attached -OH group (qualifies as a long-chain fatty alcohol)."
    
    # If no suitable chain was found, report the longest chain length encountered.
    if max_path_length_found < 13:
        return False, f"No carbon chain of at least 13 carbons found (longest found was {max_path_length_found})."
    else:
        return False, f"Found carbon chains up to {max_path_length_found} atoms, but none with an attached -OH group in a chain of length 13-22."

# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = [
        "OCCCCCCCCCC/C=C/CCCCCCCC",  # 11E-Eicosen-1-ol: 20 carbons
        "CCCCCCCCCCCCCCCCCCCCCO",    # docosan-1-ol: 22 carbons
        "CC(O)CCCCCCCCCC",            # too short (maybe 11 carbons)
    ]
    for smi in test_smiles:
        result, msg = is_long_chain_fatty_alcohol(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {msg}\n")