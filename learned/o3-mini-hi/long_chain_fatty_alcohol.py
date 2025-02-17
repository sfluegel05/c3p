"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: Long-chain fatty alcohol
Definition: A fatty alcohol with a carbon chain length ranging from C13 to C22.
A valid long-chain fatty alcohol should have a contiguous non-cyclic (acyclic) alkyl chain
of carbon atoms numbering between 13 and 22, and at least one of the carbons in that chain
must bear an –OH substituent.
Note: This is a heuristic method.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    The heuristic looks for a contiguous acyclic carbon chain (i.e. outside of rings) of length
    between 13 and 22 that has at least one carbon with an attached –OH group.
    Also, the overall molecular weight should be in the typical range for fatty alcohols.

    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule qualifies as a long-chain fatty alcohol, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight (most fatty alcohols have MW < ~500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too high for a simple long-chain fatty alcohol"
    
    # Add explicit hydrogens to help locate –OH groups.
    molH = Chem.AddHs(mol)
    
    # Build a dictionary for candidate alkyl chain carbons:
    # Only include carbons that are not part of any ring.
    # Also note for each such carbon whether it has a directly attached –OH group.
    alkyl_carbons = {}
    for atom in molH.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.IsInRing():
            # Check if any neighbor is an oxygen that is part of a hydroxyl (-OH)
            # (i.e. oxygen bonded to at least one hydrogen)
            attached_OH = False
            for nbr in atom.GetNeighbors():
                # Look for oxygen atoms
                if nbr.GetAtomicNum() == 8:
                    # If oxygen has at least one hydrogen attached, count as alcohol O.
                    if any(n.GetAtomicNum() == 1 for n in nbr.GetNeighbors()):
                        attached_OH = True
                        break
            alkyl_carbons[atom.GetIdx()] = attached_OH

    if not alkyl_carbons:
        return False, "No non-ring carbon atoms found to form an alkyl chain."
    
    # Build the connectivity graph among these acyclic (alkyl) carbons.
    # The graph is a dictionary mapping atom index -> list of neighboring atom indices.
    carbon_graph = {idx: [] for idx in alkyl_carbons}
    for bond in molH.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            idx1 = a1.GetIdx()
            idx2 = a2.GetIdx()
            if idx1 in alkyl_carbons and idx2 in alkyl_carbons:
                carbon_graph[idx1].append(idx2)
                carbon_graph[idx2].append(idx1)
    
    # Partition the graph into connected components.
    components = []
    visited_global = set()
    for node in carbon_graph.keys():
        if node not in visited_global:
            # simple DFS to get all nodes of this component.
            comp = set()
            stack = [node]
            while stack:
                current = stack.pop()
                if current in comp:
                    continue
                comp.add(current)
                for nbr in carbon_graph.get(current, []):
                    if nbr not in comp:
                        stack.append(nbr)
            components.append(list(comp))
            visited_global |= comp

    # Helper: In a tree component the (unique) simple path between two vertices can be found via DFS.
    def find_path(comp_nodes, graph, start, goal):
        # returns list of nodes representing the path (including both ends), or None if not found.
        stack = [(start, [start])]
        visited_local = set()
        while stack:
            (current, path) = stack.pop()
            if current == goal:
                return path
            visited_local.add(current)
            for nbr in graph.get(current, []):
                if nbr in comp_nodes and nbr not in path:
                    new_path = path + [nbr]
                    stack.append((nbr, new_path))
        return None

    # Now, search for any simple path within any connected component of alkyl carbons that:
    #   - Has a length (number of carbons) between 13 and 22 (inclusive)
    #   - Contains at least one carbon that directly bears a –OH group.
    candidate_found = False
    candidate_path_length = 0
    candidate_comp_info = ""
    max_path_length_overall = 0

    for comp in components:
        # Only consider components that have at least one potential candidate carbon.
        if len(comp) < 13:
            if len(comp) > max_path_length_overall:
                max_path_length_overall = len(comp)
            continue
        # For all pairs in the component (since our graph is a tree or forest, the path is unique)
        for i, start in enumerate(comp):
            for goal in comp[i+1:]:
                path = find_path(set(comp), carbon_graph, start, goal)
                if path is None:
                    continue
                path_len = len(path)
                if path_len > max_path_length_overall:
                    max_path_length_overall = path_len
                if 13 <= path_len <= 22:
                    # Check that at least one carbon in the path has an -OH substituent.
                    if any(alkyl_carbons.get(atom_idx, False) for atom_idx in path):
                        candidate_found = True
                        candidate_path_length = path_len
                        candidate_comp_info = (
                            f"Found an acyclic carbon chain of length {path_len} containing an -OH substituent."
                        )
                        break
            if candidate_found:
                break
        if candidate_found:
            break

    if candidate_found:
        return True, candidate_comp_info
    else:
        # No valid chain found. Report the maximum chain length seen.
        if max_path_length_overall < 13:
            return False, f"No acyclic carbon chain of at least 13 carbons found (longest was {max_path_length_overall})."
        else:
            return False, f"Found acyclic carbon chains up to {max_path_length_overall} atoms, but none with an -OH group on a chain between 13 and 22 carbons."

# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = [
        "OCCCCCCCCCC/C=C/CCCCCCCC",  # 11E-Eicosen-1-ol: should qualify
        "CCCCCCCCCCCCCCCCCCCCCO",    # docosan-1-ol: 22 carbons, qualifies if terminal -OH
        "CC(O)CCCCCCCCCC",            # Likely too short (~11 carbons)
        "O1[C@@H]2[C@H](O)C(/C=C/[C@H](O)CCCCCCCC)=C([C@H]([C@H]12)O)CO",  # Phomopoxide D (should qualify)
    ]
    for smi in test_smiles:
        result, msg = is_long_chain_fatty_alcohol(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {msg}\n")