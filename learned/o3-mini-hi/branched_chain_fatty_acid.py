"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
"""
Classifies: branched-chain fatty acid (BCFA)
Definition:
  A branched–chain fatty acid is a fatty acid (i.e. a carboxylic acid with an alkyl chain)
  in which the parent hydrocarbon chain (the fatty acyl chain) carries one or more extra 
  alkyl substituents. Such substituents are most commonly small (typically a methyl group),
  though occasionally slightly larger branches occur.
  
This implementation:
  1. Searches for a carboxylic acid group using two SMARTS patterns (for the protonated and 
     deprotonated forms).
  2. Identifies the fatty acyl chain as the longest acyclic chain (over carbons only) extending 
     from the acid carbon.
  3. Requires that this chain be dominant (i.e. it covers a significant fraction of the molecule’s 
     carbons) and reasonably long.
  4. Scans the main chain for a branch: for each carbon in the main chain (excluding the acid carbon),
     it checks for an extra carbon substituent. However, instead of accepting any branch it uses a DFS 
     restricted to carbons outside the chosen main chain to count the number of connected atoms. Only 
     very short branches (one or two carbons) are allowed.
     
Note: This heuristic is still approximate.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

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
    
    # Step 1. Look for a carboxylic acid group.
    # Try both protonated and deprotonated forms.
    acid_smarts_list = ["[CX3](=O)[OX2H]", "[CX3](=O)[OX2-]"]
    acid_matches = []
    for smt in acid_smarts_list:
        patt = Chem.MolFromSmarts(smt)
        acid_matches = mol.GetSubstructMatches(patt)
        if acid_matches:
            break
    if not acid_matches:
        return False, "No carboxylic acid group found; not a fatty acid."
    
    # Use the first match. In our SMARTS the first atom is the acid carbon.
    acid_idx = acid_matches[0][0]
    acid_atom = mol.GetAtomWithIdx(acid_idx)
    
    # Step 2. Identify the fatty acyl chain:
    # Get the R group attached to the acid carbon (ignoring oxygen neighbors).
    rgroup_atom = None
    for nbr in acid_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6:  # a carbon neighbor
            rgroup_atom = nbr
            break
    if rgroup_atom is None:
        return False, "No alkyl chain (R group) attached to the carboxyl carbon; not a fatty acid."
    
    # Get all carbons in the molecule.
    all_carbons = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(all_carbons) < 3:
        return False, "Too few carbon atoms to be considered a fatty acid."
    
    # Build a simple carbon-only graph
    carbon_graph = {}
    for idx in all_carbons:
        atom = mol.GetAtomWithIdx(idx)
        carbon_graph[idx] = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    
    # DFS helper to find the longest acyclic chain starting from a given atom.
    def dfs_longest(start, visited):
        best_path = [start]
        for nbr in carbon_graph.get(start, []):
            if nbr in visited:
                continue
            candidate = dfs_longest(nbr, visited | {nbr})
            if len(candidate) + 1 > len(best_path):
                best_path = [start] + candidate
        return best_path

    # Start DFS from the R-group carbon; exclude the acid carbon from the chain.
    start_idx = rgroup_atom.GetIdx()
    main_chain = dfs_longest(start_idx, {acid_idx, start_idx})
    
    # The complete fatty acyl unit includes the acid carbon.
    acyl_chain_set = set(main_chain) | {acid_idx}
    
    # Step 3a. Ensure the main chain is long enough.
    if len(main_chain) < 1:
        return False, "The fatty acyl chain appears too short."
    
    # Step 3b. Check that the acyl chain is the dominant carbon skeleton.
    ratio = len(acyl_chain_set) / len(all_carbons)
    if len(all_carbons) > 8 and ratio < 0.5:
        return False, "The fatty acyl chain is not the dominant carbon skeleton in the molecule."
    
    # Step 4. Check for presence of a small branch on the main chain (excludes the acid carbon).
    # Here a branch must be small (size 1 or 2 carbons) to be accepted.
    max_branch_atoms = 2  # allow methyl (1) or ethyl (2); larger branches cause false positives
    
    # helper: given a branch starting at a carbon (which is not in the main chain),
    # traverse all connected carbons (restricted to those not in main chain) and count them.
    def get_branch_size(start, exclude_set):
        stack = [start]
        branch_atoms = set()
        while stack:
            current = stack.pop()
            if current in branch_atoms:
                continue
            branch_atoms.add(current)
            atom = mol.GetAtomWithIdx(current)
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx in exclude_set or nbr_idx in branch_atoms:
                        continue
                    stack.append(nbr_idx)
        return len(branch_atoms)
    
    branch_found = False
    branch_reason = ""
    # For each carbon in the main chain (we skip acid carbon because it is considered the reactive head)
    for atom_idx in main_chain:
        atom = mol.GetAtomWithIdx(atom_idx)
        # Look for neighboring carbons that are not part of the fatty acyl chain.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            nbr_idx = nbr.GetIdx()
            if nbr_idx in acyl_chain_set:
                continue
            # Found a candidate branch; compute its size.
            branch_size = get_branch_size(nbr_idx, acyl_chain_set)
            # Accept if branch is small.
            if branch_size <= max_branch_atoms:
                branch_found = True
                branch_reason = f"Found a branch of {branch_size} carbon(s) off the main chain."
                break
        if branch_found:
            break

    if not branch_found:
        return False, "No small alkyl substituents (branch) found on the fatty acyl chain."
    
    return True, "Contains a carboxylic acid group with a dominant fatty acyl chain and a small alkyl branch. " + branch_reason


# Example usage (for quick testing):
if __name__ == "__main__":
    examples = [
        # True positives:
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
        ("CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\C(O)=O", "(2E,6E,10E)-geranylgeranic acid"),
        ("C[C@@H](O)CC[C@@H](CC(O)=O)C(C)=C", "(3S,6R)-6-hydroxy-3-isopropenylheptanoic acid"),
        ("CCCCCCCCCC(C)CC(O)=O", "3-methylundecanoic acid"),
        ("OC(=O)C(CC=C)(C)C", "2,2-dimethyl-4-pentenoic acid"),
        ("OC(=O)CC(CC)(C)C", "beta,beta-dimethyl valeric acid"),
        ("OC(=O)C(CCC)CC", "alpha-ethyl valeric acid"),
        ("CCCCCC[C@H](C)[C@@H](O)CC(O)=O", "(3S,4S)-3-hydroxy-4-methyldecanoic acid"),
        ("CC(C)C(O)=O", "isobutyric acid"),
        ("CC(C)CCCCCCCCCCCCCCCCCCCCCCCCC(O)=O", "26-methylheptacosanoic acid"),
        ("OC(=O)CC(C)C=C", "3-methyl-4-pentenoic acid"),
        ("O[C@H]([C@@H](CC)C)C(O)=O", "(2R,3R)-2-hydroxy-3-methylpentanoic acid"),
        ("CC(C)C[C@H](O)C(O)=O", "(S)-2-hydroxy-4-methylpentanoic acid"),
        ("CCCCCCCCCCCCCCCCCCCCCCCC[C@H]([C@H](O)CCCCCCCCCCC[C@H]1C[C@H]1CCCCCCCCCCCCCC[C@H]1C[C@H]1CCCCCCCCCCCCCCCCCCCC)C(O)=O",
         "(2R)-2-{(1R)-1-hydroxy-12-[(1S,2R)-2-{14-[(1S,2R)-2-icosylcyclopropyl]tetradecyl}cyclopropyl]dodecyl}hexacosanoic acid"),
        # False positive candidates:
        ("OC(=O)CCC(CCCC)CC", "4-Ethyloctanoic acid"),
        # False negatives:
        ("[H][C@@]12[C@H](CCN1CC=C2COC(=O)[C@](O)([C@H](C)O)C(C)(C)O)OC(=O)C(\\C)=C/C", "heliosupine"),
        ("[O-]C(=O)CCC([N+](C)(C)C)C", "4-aminovaleric acid betaine")
    ]
    
    for smi, name in examples:
        result, reason = is_branched_chain_fatty_acid(smi)
        print(f"SMILES: {smi}\nNAME: {name}\nResult: {result}\nReason: {reason}\n{'-'*70}")