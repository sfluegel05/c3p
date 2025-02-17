"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: Long-chain fatty alcohol
Definition: A fatty alcohol with a contiguous carbon chain of length C13–C22 (relaxed to C12 if rings are present)
that carries at least one hydroxyl (–OH) group (not part of a carboxylic acid) and does not contain extra acid functionality.
Additional heuristics:
  – Overall molecular weight must be below ~500 Da.
  – The candidate chain should account for most of the molecule’s carbon atoms (~80%).
Note: This heuristic method is approximate.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule qualifies as a long-chain fatty alcohol.
    It demands that (a) the molecular weight is below 500 Da,
    (b) the molecule does not have carboxylic acid/carboxylate groups,
    (c) there is a long contiguous carbon chain (as measured by the longest simple path in the carbon network)
        of length between 13 and 22 (or lower bound of 12 if rings are present),
    (d) that chain covers at least 80% of the total carbons,
    and (e) at least one carbon on the candidate chain bears a valid hydroxyl group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule qualifies as a long-chain fatty alcohol, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules with a high molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too high for a simple long-chain fatty alcohol"
    
    # Work on a version with explicit hydrogens.
    molH = Chem.AddHs(mol)
    
    # Reject if the molecule contains a carboxylic acid group.
    # SMARTS for acid: either protonated [CX3](=O)[OX2H] or deprotonated [CX3](=O)[O-]
    acid_smarts = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    acid_smarts2 = Chem.MolFromSmarts("[CX3](=O)[O-]")
    if (acid_smarts and molH.HasSubstructMatch(acid_smarts)) or (acid_smarts2 and molH.HasSubstructMatch(acid_smarts2)):
        return False, "Contains a carboxylic acid group (or carboxylate) which is not allowed in a simple fatty alcohol"
    
    # Build the complete carbon graph (all C atoms, all C–C bonds regardless of ring membership)
    carbon_nodes = {}
    for atom in molH.GetAtoms():
        if atom.GetAtomicNum() == 6:
            carbon_nodes[atom.GetIdx()] = atom
    if not carbon_nodes:
        return False, "No carbon atoms for chain analysis."
        
    # Build graph: For every bond between two carbon atoms, add an undirected edge.
    carbon_graph = {idx: [] for idx in carbon_nodes}
    for bond in molH.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            i1 = a1.GetIdx()
            i2 = a2.GetIdx()
            if i1 in carbon_graph and i2 in carbon_graph:
                carbon_graph[i1].append(i2)
                carbon_graph[i2].append(i1)
    
    total_carbons = len(carbon_nodes)
    
    # Use DFS to get the longest simple (nonrepeating) path among carbon atoms.
    # For our typically small molecules this brute force search is acceptable.
    best_path = []
    # Use recursion with backtracking.
    def dfs(current, path, visited):
        nonlocal best_path
        # update best_path if current path is longer
        if len(path) > len(best_path):
            best_path = path[:]
        for nbr in carbon_graph.get(current, []):
            if nbr not in visited:
                visited.add(nbr)
                path.append(nbr)
                dfs(nbr, path, visited)
                path.pop()
                visited.remove(nbr)
    
    for start in carbon_nodes:
        dfs(start, [start], {start})
    
    candidate_length = len(best_path)
    
    # Determine the lower bound on chain length.
    # When rings are present, the longest simple path using all carbons may be slightly shorter,
    # so we allow a lower bound of 12 instead of 13.
    if molH.GetRingInfo().NumRings() > 0:
        lower_bound = 12
    else:
        lower_bound = 13
        
    if candidate_length < lower_bound:
        return False, f"No contiguous carbon chain of at least {lower_bound} carbons found (longest was {candidate_length})."
    if candidate_length > 22:
        return False, f"Longest carbon chain (length {candidate_length}) exceeds the desired upper limit of 22."
    
    # Dominance criterion: the candidate chain should cover at least 80% of the molecule's carbon atoms.
    if candidate_length < 0.8 * total_carbons:
        return False, (f"Longest chain length ({candidate_length}) accounts for less than 80% of all {total_carbons} carbon atoms, "
                       "indicating a complex framework.")
    
    # Check that at least one carbon in the candidate chain carries a valid hydroxyl (–OH) group.
    # We look for an oxygen (atomic num 8) that has at least one hydrogen and is attached to a carbon in the candidate.
    has_valid_OH = False
    for idx in best_path:
        carbon = carbon_nodes[idx]
        for nbr in carbon.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                # require at least one hydrogen attached to the oxygen
                if any(n.GetAtomicNum() == 1 for n in nbr.GetNeighbors()):
                    has_valid_OH = True
                    break
        if has_valid_OH:
            break
    if not has_valid_OH:
        return False, ("Candidate carbon chain (length {}) does not have any valid hydroxyl (-OH) group attached."
                       .format(candidate_length))
    
    msg = (f"Found a contiguous carbon chain of length {candidate_length} "
           f"(covering {total_carbons} total carbons) with a valid -OH substituent.")
    return True, msg

# Example usage (for testing)
if __name__ == "__main__":
    test_smiles = [
        "OCCCCCCCCCC/C=C/CCCCCCCC",  # 11E-Eicosen-1-ol; should qualify.
        "CCCCCCCCCCCCCCCCCCCCCO",    # docosan-1-ol; 22-carbons chain with terminal -OH.
        "CC(O)CCCCCCCCCC",            # Likely too short.
        "O1[C@@H]2[C@H](O)C(/C=C/[C@H](O)CCCCCCCC)=C([C@H]([C@H]12)O)CO",  # Phomopoxide D; candidate chain might be 12 => allowed if rings present.
    ]
    for smi in test_smiles:
        res, reason = is_long_chain_fatty_alcohol(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n")