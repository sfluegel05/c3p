"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
"""
Classifies: branched-chain fatty acid (BCFA)
Definition:
  A branched–chain fatty acid is a fatty acid (i.e. a carboxylic acid with an alkyl chain)
  in which the parent hydrocarbon chain (the fatty acyl chain) carries one or more extra 
  alkyl substituents. Such substituents are commonly small (a methyl group), although
  others may occur.
  
This implementation:
  1. Looks for a carboxylic acid moiety (protonated or deprotonated).
  2. Identifies the R group (a carbon neighbor of the acid carbon).
  3. Uses depth–first search (DFS) over strictly carbon–only, acyclic paths (to avoid cycles)
     to find the longest chain (the main fatty acyl backbone).
  4. Checks that this chain is:
       a. Long enough (or, if the molecule is very small, is acceptable),
       b. Dominant among the molecule’s carbon atoms (if the molecule is large), and 
       c. Contains at least one branch: meaning that one of the atoms in the main chain 
          carries a carbon substituent that is not part of the main chain (nor the acid C).
  
Note: This heuristic may not be perfect. Improvements in branch‐identification or handling 
      more complex molecules would require more advanced graph algorithms.
"""
from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid (BCFA) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule qualifies as a BCFA, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Step 1. Find a carboxylic acid group.
    # The SMARTS here accepts a carbonyl carbon attached to an -OH or -O(-).
    acid_smarts = Chem.MolFromSmarts("[CX3](=O)[O;H1,O-]")
    acid_matches = mol.GetSubstructMatches(acid_smarts)
    if not acid_matches:
        return False, "No carboxylic acid group found; not a fatty acid."
    
    # Assume the first match corresponds to the fatty acid head.
    # In the SMARTS the first atom is the acid carbon.
    acid_idx = acid_matches[0][0]
    acid_atom = mol.GetAtomWithIdx(acid_idx)
    
    # Step 2. Identify the R group (the alkyl chain attached to the acid carbon).
    # We ignore oxygen neighbors (which are part of the acid group).
    rgroup_atom = None
    for nbr in acid_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6:  # carbon neighbor
            rgroup_atom = nbr
            break
    if rgroup_atom is None:
        return False, "No alkyl chain (R group) attached to the carboxyl carbon; not a fatty acid."

    # Create a list of all carbon atom indices in the molecule.
    all_carbons = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(all_carbons) < 3:
        return False, "Too few carbon atoms to be considered a fatty acid."
    
    # Step 3. Build a local subgraph over carbons.
    carbon_graph = {}
    for idx in all_carbons:
        atom = mol.GetAtomWithIdx(idx)
        # Consider only neighboring carbons.
        carbon_graph[idx] = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    
    # DFS helper to find the longest acyclic chain (list of carbon indices) from a starting atom.
    def dfs_longest(start, visited):
        best_path = [start]
        for nbr in carbon_graph.get(start, []):
            if nbr in visited:
                continue
            candidate = dfs_longest(nbr, visited | {nbr})
            if len(candidate) + 1 > len(best_path):
                best_path = [start] + candidate
        return best_path

    # Start DFS from the rgroup_atom.
    start_idx = rgroup_atom.GetIdx()
    main_chain = dfs_longest(start_idx, {acid_idx, start_idx})  # Avoid going back to acid carbon.
    
    # Include the acid carbon in the acyl unit for ratio calculations.
    acyl_chain_set = set(main_chain) | {acid_idx}
    
    # Step 4a. Check that the chain is long enough.
    # We allow very short chains if the overall molecule is small.
    if len(main_chain) < 1:
        return False, "The fatty acyl chain appears too short."
    
    # Step 4b. Check that the acyl chain accounts for a significant part of the molecule.
    ratio = len(acyl_chain_set) / len(all_carbons)
    # For small molecules (total carbons <= 8) we relax this threshold.
    if len(all_carbons) > 8 and ratio < 0.5:
        return False, "The fatty acyl chain is not the dominant carbon skeleton in the molecule."
    
    # Step 5. Check for the presence of at least one branch.
    # For each carbon atom in the main chain, check its carbon neighbors.
    branch_found = False
    for atom_idx in main_chain:
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            nbr_idx = nbr.GetIdx()
            # Skip if neighbor is in the main chain or is the acid carbon.
            if nbr_idx in acyl_chain_set:
                continue
            # If found an extra carbon substituent, we mark as branch.
            branch_found = True
            break
        if branch_found:
            break

    if not branch_found:
        return False, "No alkyl substituents (branch) found on the fatty acyl chain."
    
    return True, "Contains a carboxylic acid group with a fatty acyl chain that is dominant and carries an alkyl branch."

# Example usage (for quick testing):
if __name__ == "__main__":
    test_examples = [
        # True positives (known BCFAs)
        ("CC\\C(CC[C@H]1O[C@@]1(C)CC)=C/CC\\C(C)=C\\C(O)=O", "juvenile hormone I acid"),
        ("CC(C)=CCCC(C)=CCCC(C)=CC(O)=O", "farnesoic acid"),
        ("CCCCCCCCCCCCCCCCCCC(C)CC(C)\\C=C(/C)C(O)=O", "(E)-2,4,6-trimethyltetracos-2-enoic acid"),
        ("CC(C)CCCCCCCCCCCCCCCCCCC(O)=O", "20-methylhenicosanoic acid"),
        ("CCCC(C)C(O)=O", "2-methylvaleric acid"),
        ("CC(CCCCCCC/C=C/C(=O)O)C", "(E)-11-methyldodec-2-enoic acid"),
        ("OC(=O)C(CC)(CC)C", "2-ethyl-2-methyl-butanoic acid"),
        ("CC(C)C(=O)CC(O)=O", "4-methyl-3-oxopentanoic acid"),
        ("CC(C)C[C@@H](O)C(O)=O", "(R)-2-hydroxy-4-methylpentanoic acid"),
        ("OC(=O)/C=C(\\CC)/C", "3-methyl-2Z-pentenoic acid"),
        ("C(C(C(O)=O)O)C(C)C", "2-hydroxy-4-methylvaleric acid"),
        ("CCC\\C(=C/CC)C(O)=O", "2-n-Propyl-2-pentenoic acid"),
        # A known false positive candidate (should return False)
        ("OC(=O)CCC(CCCC)CC", "4-Ethyloctanoic acid"),
        # A known false negative example (heavily functionalized molecule)
        ("[H][C@@]12[C@H](CCN1CC=C2COC(=O)[C@](O)([C@H](C)O)C(C)(C)O)OC(=O)C(\\C)=C/C", "heliosupine"),
        # isobutyric acid (should be accepted)
        ("CC(C)C(O)=O", "isobutyric acid")
    ]
    
    for smi, name in test_examples:
        result, reason = is_branched_chain_fatty_acid(smi)
        print(f"SMILES: {smi}\nNAME: {name}\nResult: {result}\nReason: {reason}\n{'-'*60}")