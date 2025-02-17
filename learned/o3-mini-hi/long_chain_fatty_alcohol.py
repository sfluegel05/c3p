"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: Long-chain fatty alcohol
Definition: A fatty alcohol with a carbon chain length ranging from C13 to C22.
A valid long-chain fatty alcohol should have a contiguous acyclic (non‐ring) alkyl chain
of carbon atoms numbering between 13 and 22, and at least one of the carbons in that chain
must bear a –OH substituent (and that –OH should not simply be part of a carboxylic acid).
Also, the overall molecule should be “simple” (e.g. MW < ~500 Da, and the candidate chain 
comprises most of the available non‐ring carbons).
Note: This is a heuristic method.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule qualifies as a long-chain fatty alcohol based on its SMILES string.
    The heuristic looks for a contiguous acyclic carbon chain (i.e. outside of rings) of length between 13 and 22
    that shows at least one hydroxyl (–OH) substituent (not simply from a carboxylic acid).
    Also, the molecule’s overall molecular weight should be below ~500 Da and the candidate chain should majorly
    account for the non‐ring carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule qualifies as a long-chain fatty alcohol, False otherwise.
        str: Explanation or reason for the classification.
    """
    # Parse SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check overall molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too high for a simple long-chain fatty alcohol"
    
    # We work with an explicit hydrogen copy to help judge –OH groups.
    molH = Chem.AddHs(mol)
    
    # Helper: Check if a given oxygen (in molH) is part of a carboxylic acid:
    # Look for an oxygen that is attached to a carbon that is double bonded to an oxygen.
    def is_acidic_O(oxygen):
        # oxygen must have at least one hydrogen (for -OH) and be bonded to a carbon that is also double-bonded to an O.
        if oxygen.GetAtomicNum() != 8:
            return False
        if not any(neigh.GetAtomicNum() == 1 for neigh in oxygen.GetNeighbors()):
            return False
        for nbr in oxygen.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                for bond in nbr.GetBonds():
                    # if the bond is double and the other end is oxygen then this oxygen is likely part of a carboxylic acid
                    if bond.GetBondTypeAsDouble() == 2.0:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            return True
        return False

    # Build dictionary of individual non‐ring carbon atoms.
    # For each such carbon mark if it has at least one oxygen neighbor that is a hydroxyl moiety.
    alkyl_carbons = {}
    for atom in molH.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.IsInRing():
            has_alcoh_OH = False
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    # Check that this oxygen is –OH (has a hydrogen) and is not part of a carboxylic acid.
                    if any(n.GetAtomicNum() == 1 for n in nbr.GetNeighbors()):
                        if not is_acidic_O(nbr):
                            has_alcoh_OH = True
                            break
            alkyl_carbons[atom.GetIdx()] = has_alcoh_OH

    if not alkyl_carbons:
        return False, "No non‐ring carbon atoms found to form an alkyl chain."
    
    # Build a connectivity graph among these acyclic carbons.
    # Graph: atom_idx -> list of neighbor atom indices (only for non‐ring carbons)
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
    visited_overall = set()
    for node in carbon_graph:
        if node not in visited_overall:
            comp = set()
            stack = [node]
            while stack:
                curr = stack.pop()
                if curr in comp:
                    continue
                comp.add(curr)
                for nbr in carbon_graph[curr]:
                    if nbr not in comp:
                        stack.append(nbr)
            components.append(list(comp))
            visited_overall |= comp

    # Helper: Find the longest simple path in a component.
    # (Since components are small, we use an exhaustive DFS over all pairs.)
    def longest_path_in_component(comp, graph):
        best_path = []
        comp_set = set(comp)
        # For each pair in the component, find the simple (non-repeating) path.
        # (For efficiency, we start DFS at each node.)
        def dfs(current, target, path, visited):
            nonlocal best_path
            if current == target:
                if len(path) > len(best_path):
                    best_path = path[:]
                return
            for nbr in graph.get(current, []):
                if nbr in comp_set and nbr not in visited:
                    visited.add(nbr)
                    path.append(nbr)
                    dfs(nbr, target, path, visited)
                    path.pop()
                    visited.remove(nbr)
        comp_nodes = list(comp_set)
        # Try every pair: 
        for i in range(len(comp_nodes)):
            for j in range(i, len(comp_nodes)):
                start = comp_nodes[i]
                target = comp_nodes[j]
                dfs(start, target, [start], {start})
        return best_path

    # Look through each connected component for a candidate chain.
    candidate_found = False
    candidate_msg = ""
    # Also record total non-ring carbon count.
    total_nonring = len(alkyl_carbons)
    best_path_overall = []
    
    for comp in components:
        # Only consider if the component is long enough to possibly have a chain.
        if len(comp) < 13:
            continue
        lp = longest_path_in_component(comp, carbon_graph)
        if len(lp) > len(best_path_overall):
            best_path_overall = lp
        if 13 <= len(lp) <= 22:
            # Check if at least one carbon in this candidate path bears an -OH (nonacidic)
            if any(alkyl_carbons.get(idx, False) for idx in lp):
                # Also enforce that the candidate chain is “large” relative to the total non‐ring carbons.
                # (for a simple fatty alcohol the chain should dominate the structure)
                if len(lp) >= 0.8 * total_nonring:
                    candidate_found = True
                    candidate_msg = (f"Found an acyclic carbon chain of length {len(lp)} (out of {total_nonring} non‐ring carbons) "
                                     "with a hydroxyl substituent not part of an acid group.")
                    break
    if candidate_found:
        return True, candidate_msg
    else:
        if best_path_overall:
            if len(best_path_overall) < 13:
                return False, f"No acyclic carbon chain of at least 13 carbons found (longest was {len(best_path_overall)})."
            else:
                return False, (f"Found an acyclic carbon chain of length {len(best_path_overall)} "
                               f"but did not meet all criteria (e.g. proper -OH placement or chain dominance).")
        else:
            return False, "No candidate acyclic carbon chain found."

# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = [
        "OCCCCCCCCCC/C=C/CCCCCCCC",  # 11E-Eicosen-1-ol; should qualify
        "CCCCCCCCCCCCCCCCCCCCCO",    # docosan-1-ol: 22 carbons, qualifies if terminal -OH
        "CC(O)CCCCCCCCCC",            # Likely too short (~11 carbons)
        "O1[C@@H]2[C@H](O)C(/C=C/[C@H](O)CCCCCCCC)=C([C@H]([C@H]12)O)CO",  # Phomopoxide D; should now qualify
    ]
    for smi in test_smiles:
        result, msg = is_long_chain_fatty_alcohol(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {msg}\n")