"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: Long-chain fatty alcohol
Definition: A fatty alcohol with a carbon chain length ranging from C13 to C22.
The structure should have a contiguous (acyclic, i.e. through non‐ring bonds) 
carbon chain in that range that carries at least one hydroxyl group (–OH) 
which is not part of a carboxylic acid.
Additional heuristics:
  – Overall molecular weight should be below ~500 Da.
  – The candidate chain should account for most of the (acyclic) carbon atoms.
  – Molecules containing carboxylic acid groups are rejected.
Note: This heuristic method is approximate.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule qualifies as a long-chain fatty alcohol.
    It requires that the molecule has a contiguous chain of carbons (via non-ring bonds)
    with a length between 13 and 22, that chain includes at least one –OH group not part 
    of a carboxylic acid, the overall molecular weight is <500 Da, and the chain dominates 
    the acyclic carbon skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule qualifies as a long-chain fatty alcohol, False otherwise.
        str: Explanation or reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject if molecular weight is too high
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too high for a simple long-chain fatty alcohol"
    
    # Create an explicit-hydrogens copy to better judge –OH groups.
    molH = Chem.AddHs(mol)
    
    # Pre-filter: reject molecules containing carboxylic acid groups.
    # SMARTS for carboxylic acid: a carbon double bonded to oxygen and bonded to an -OH.
    acid_smarts = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    if acid_smarts is not None and molH.HasSubstructMatch(acid_smarts):
        return False, "Contains a carboxylic acid group which is not allowed in a simple fatty alcohol"
    
    # Helper function: Check if an oxygen (in molH) is acidic.
    def is_acidic_O(oxygen):
        # Must be oxygen with at least one bound hydrogen.
        if oxygen.GetAtomicNum() != 8:
            return False
        if not any(neigh.GetAtomicNum() == 1 for neigh in oxygen.GetNeighbors()):
            return False
        # Look at each neighbor carbon; if that carbon is double bonded to another O then mark it acidic.
        for nbr in oxygen.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                for bond in nbr.GetBonds():
                    if bond.GetBondTypeAsDouble() == 2.0:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            return True
        return False

    # Build a connectivity graph of carbon atoms (all carbons)
    # but only use bonds that are not part of a ring.
    carbon_nodes = {}  # map: atom index -> atom
    for atom in molH.GetAtoms():
        if atom.GetAtomicNum() == 6:
            carbon_nodes[atom.GetIdx()] = atom
    
    # Build graph: edges only if bond is between two carbons and bond.IsInRing() is False.
    carbon_graph = {idx: [] for idx in carbon_nodes}
    for bond in molH.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            # Only add edge if this bond is not in a ring.
            if not bond.IsInRing():
                i1 = a1.GetIdx()
                i2 = a2.GetIdx()
                if i1 in carbon_graph and i2 in carbon_graph:
                    carbon_graph[i1].append(i2)
                    carbon_graph[i2].append(i1)
    
    if not carbon_graph:
        return False, "No carbon atoms available for chain analysis."
    
    total_carbons = len(carbon_graph)  # Number of carbon atoms in the acyclic connectivity graph

    # Partition the graph into connected components.
    visited_overall = set()
    components = []
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

    # Helper: Find the longest simple path in a component by DFS (exhaustive search).
    def longest_path_in_component(comp, graph):
        best_path = []
        comp_set = set(comp)
        # Recursive DFS for simple paths.
        def dfs(current, path, visited):
            nonlocal best_path
            # update best_path if this path is longer
            if len(path) > len(best_path):
                best_path = path[:]
            for nbr in graph.get(current, []):
                if nbr in comp_set and nbr not in visited:
                    visited.add(nbr)
                    path.append(nbr)
                    dfs(nbr, path, visited)
                    path.pop()
                    visited.remove(nbr)
        for start in comp:
            dfs(start, [start], {start})
        return best_path

    candidate_found = False
    candidate_chain = []
    candidate_msg = ""
    
    # Check each connected component for a candidate chain.
    # We require that the candidate chain length must be between 13 and 22.
    # Also, we require that it accounts for at least 80% of the carbons in the graph.
    for comp in components:
        lp = longest_path_in_component(comp, carbon_graph)
        # if the longest found in this component is too long, we may consider subpaths,
        # but here we simply check if the longest path length falls in our range.
        if 13 <= len(lp) <= 22:
            # Check that at least one carbon in the candidate has a –OH (nonacidic) attached.
            has_valid_OH = False
            for idx in lp:
                carbon_atom = carbon_nodes[idx]
                for nbr in carbon_atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 8:
                        # Must have a hydrogen attached and not be acidic.
                        if any(n.GetAtomicNum() == 1 for n in nbr.GetNeighbors()) and not is_acidic_O(nbr):
                            has_valid_OH = True
                            break
                if has_valid_OH:
                    break
            if not has_valid_OH:
                # If candidate chain lacks a suitable -OH then skip this component.
                continue
            # Check that the candidate chain is “dominant”
            if len(lp) < 0.8 * total_carbons:
                # Even if we have a chain of correct length, if it does not account for most carbons,
                # the molecule is likely too complex.
                continue
            # Found candidate!
            candidate_found = True
            candidate_chain = lp
            candidate_msg = (f"Found an acyclic carbon chain of length {len(lp)} "
                             f"(out of {total_carbons} connected carbons) with a valid -OH substituent.")
            break

    if candidate_found:
        return True, candidate_msg
    else:
        # If no candidate chain met criteria, report based on the longest chain overall.
        overall_best = []
        for comp in components:
            lp = longest_path_in_component(comp, carbon_graph)
            if len(lp) > len(overall_best):
                overall_best = lp
        if overall_best:
            if len(overall_best) < 13:
                return False, f"No acyclic carbon chain of at least 13 carbons found (longest was {len(overall_best)})."
            else:
                return False, (f"Found an acyclic carbon chain of length {len(overall_best)} "
                               "but it did not meet all criteria (e.g. valid -OH placement or chain dominance).")
        else:
            return False, "No candidate acyclic carbon chain found."


# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = [
        "OCCCCCCCCCC/C=C/CCCCCCCC",  # 11E-Eicosen-1-ol; should qualify
        "CCCCCCCCCCCCCCCCCCCCCO",    # docosan-1-ol; 22-carbon chain with terminal -OH
        "CC(O)CCCCCCCCCC",            # Likely too short (11 carbons)
        "O1[C@@H]2[C@H](O)C(/C=C/[C@H](O)CCCCCCCC)=C([C@H]([C@H]12)O)CO",  # Phomopoxide D; should now be captured if the chain at the ring junction is taken into account
    ]
    for smi in test_smiles:
        res, msg = is_long_chain_fatty_alcohol(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {msg}\n")