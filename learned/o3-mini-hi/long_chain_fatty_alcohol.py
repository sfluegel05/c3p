"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: Long-chain fatty alcohol
Definition: A fatty alcohol with a carbon chain length ranging from C13 to C22.
A valid long-chain fatty alcohol should have an acyclic alkyl chain containing 13–22 carbon atoms,
with at least one carbon in that chain bearing an –OH substituent.
Note: This is a heuristic method.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    We require that the molecule has a contiguous (acyclic) chain of carbon atoms of length 13 to 22,
    with at least one carbon attached to an -OH group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule qualifies as a long-chain fatty alcohol, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add hydrogens so we can detect -OH groups.
    molH = Chem.AddHs(mol)
    
    # Build a dictionary for carbon atoms:
    # Only consider acyclic carbons (not in any ring) because we want an alkyl chain.
    carbon_atoms = {}
    for atom in molH.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.IsInRing():
            # Initially mark as not having an -OH attached.
            carbon_atoms[atom.GetIdx()] = False
    
    # Determine which carbon atoms (in the acyclic set) carry an -OH substituent.
    # We look for an oxygen attached to the carbon that is bonded to exactly one hydrogen.
    for atom in molH.GetAtoms():
        if atom.GetAtomicNum() == 6 and (atom.GetIdx() in carbon_atoms):
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    # Check that this oxygen is part of an alcohol (-OH).
                    # Expect oxygen to have exactly 2 neighbors: one hydrogen and one carbon.
                    if len(nbr.GetNeighbors()) == 2:
                        hasH = any(n.GetAtomicNum() == 1 for n in nbr.GetNeighbors())
                        if hasH:
                            carbon_atoms[atom.GetIdx()] = True
                            break
    
    # Construct a graph (as a dictionary) with only the acyclic carbon atoms.
    carbon_graph = {idx: [] for idx in carbon_atoms.keys()}
    for bond in molH.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Consider only bonds between carbon atoms that are acyclic.
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            if (not a1.IsInRing()) and (not a2.IsInRing()):
                i1 = a1.GetIdx()
                i2 = a2.GetIdx()
                if i1 in carbon_graph and i2 in carbon_graph:
                    carbon_graph[i1].append(i2)
                    carbon_graph[i2].append(i1)
    
    # We now search for any simple (acyclic) path in the carbon_graph that has length in [13,22]
    # and for which at least one carbon in the path has an attached -OH group.
    found = False
    found_path_length = 0
    reason = ""
    max_path_length_found = 0  # track longest acyclic chain encountered
    
    # DFS to search for suitable acyclic carbon chain.
    def dfs(current, visited, path):
        nonlocal found, found_path_length, max_path_length_found
        path.append(current)
        visited.add(current)
        if len(path) > max_path_length_found:
            max_path_length_found = len(path)
        # Check if the current simple path is within desired range.
        if 13 <= len(path) <= 22:
            # Check if any carbon in this path has an -OH substituent.
            if any(carbon_atoms[atom_idx] for atom_idx in path):
                found = True
                found_path_length = len(path)
                visited.remove(current)
                path.pop()
                return
        # If path is already the maximum allowed, do not extend further.
        if len(path) >= 22:
            visited.remove(current)
            path.pop()
            return
        # Continue DFS for unvisited neighbors.
        for neighbor in carbon_graph.get(current, []):
            if neighbor not in visited:
                dfs(neighbor, visited, path)
                if found:
                    return
        visited.remove(current)
        path.pop()
    
    # Start DFS from each acyclic carbon atom.
    for node in list(carbon_graph.keys()):
        if found:
            break
        dfs(node, set(), [])
    
    if found:
        return True, f"Found an acyclic carbon chain of length {found_path_length} with an attached -OH group (qualifies as a long-chain fatty alcohol)."
    
    # No valid chain was found.
    if max_path_length_found < 13:
        return False, f"No acyclic carbon chain of at least 13 carbons found (longest was {max_path_length_found})."
    else:
        return False, f"Found acyclic carbon chains up to {max_path_length_found} atoms, but none with an -OH group on a chain between 13 and 22 carbons."
    
# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = [
        "OCCCCCCCCCC/C=C/CCCCCCCC",  # 11E-Eicosen-1-ol: should qualify (~20 carbons, with -OH)
        "CCCCCCCCCCCCCCCCCCCCCO",    # docosan-1-ol: 22 carbons, qualifies if acyclic and -OH at terminus
        "CC(O)CCCCCCCCCC",            # Too short (likely 11 carbons)
        # A complex aromatic molecule (should now be disregarded as it is not acyclic),
        "c1cc2cc3cc4ccc5cc6cc7ccc8cc9ccc%10cc%11cc%12ccc%13cc%14cc%15ccc%16cc1c1c2c2c3c3c4c5c4c6c5c7c8c6c9c%10c7c%11c8c%12c%13c9c%14c%10c%15c%16c1c1c2c2c3c4c3c5c6c7c4c8c9c(c%101)c2c34"
    ]
    for smi in test_smiles:
        result, msg = is_long_chain_fatty_alcohol(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {msg}\n")