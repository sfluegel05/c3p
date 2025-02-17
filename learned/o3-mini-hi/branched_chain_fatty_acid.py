"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
"""
Classifies: branched-chain fatty acid (BCFA)
Definition:
  A branched–chain fatty acid is a fatty acid (i.e. a carboxylic acid with an alkyl chain)
  in which the parent hydrocarbon chain (the fatty acyl chain) carries one or more extra 
  alkyl substituents. Such substituents are most commonly small (typically a methyl group),
  though sometimes ethyl substituents are observed.
  
This implementation:
  1. Searches for a carboxylic acid group using SMARTS patterns for both protonated and deprotonated forms.
  2. Builds a carbon–only graph of the molecule and finds the longest simple (acyclic) path that starts at the acid carbon.
     This chain is taken as the fatty acyl chain.
  3. Checks that this chain is long enough and dominates the molecule’s carbon skeleton.
  4. Scans the internal carbons of the main chain for a substituent branch (outside the chain) whose DFS‐computed size (number of connected carbons)
     is 1 or 2. Only such small branches are accepted.
     
Note: This heuristic is approximate.
"""
from rdkit import Chem
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
    
    # Step 1. Look for a carboxylic acid group using SMARTS.
    acid_smarts_list = ["[CX3](=O)[OX2H]", "[CX3](=O)[OX2-]"]
    acid_matches = []
    for smt in acid_smarts_list:
        patt = Chem.MolFromSmarts(smt)
        matches = mol.GetSubstructMatches(patt)
        if matches:
            acid_matches = matches
            break
    if not acid_matches:
        return False, "No carboxylic acid group found; not a fatty acid."
        
    # Use the first match – in our SMARTS the 0th atom is the acid carbon.
    acid_idx = acid_matches[0][0]
    
    # Step 2. Build a carbon-only graph from the molecule.
    all_carbons = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(all_carbons) < 3:
        return False, "Too few carbon atoms to be considered a fatty acid."
    
    # Build a dictionary mapping each carbon index to its carbon neighbours.
    carbon_graph = { idx: [] for idx in all_carbons }
    for idx in all_carbons:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                carbon_graph[idx].append(nbr.GetIdx())
    
    # Use DFS to enumerate all simple (acyclic) paths starting at the acid carbon.
    # We then select the longest path as the fatty acyl chain.
    def dfs_paths(current, path):
        paths = []
        extended = False
        for nbr in carbon_graph.get(current, []):
            if nbr not in path:
                extended = True
                paths.extend(dfs_paths(nbr, path + [nbr]))
        if not extended:
            return [path]
        return paths
    
    all_paths = dfs_paths(acid_idx, [acid_idx])
    if not all_paths:
        return False, "Could not determine a fatty acyl chain from the acid carbon."
    
    # Choose the longest path (if multiple paths have the same length, one is arbitrarily selected).
    main_chain = max(all_paths, key=lambda p: len(p))
    acyl_chain_set = set(main_chain)
    
    # Step 3a. Ensure the main chain is long enough (we expect at least 3 carbons in the acyl chain).
    if len(main_chain) < 3:
        return False, "The fatty acyl chain appears too short."
    
    # Step 3b. Verify that the acyl chain constitutes a dominant fraction of all carbons.
    ratio = len(acyl_chain_set) / len(all_carbons)
    if len(all_carbons) > 8 and ratio < 0.5:
        return False, "The fatty acyl chain is not the dominant carbon skeleton in the molecule."
    
    # Step 4. Look for a small branch (alkyl substituent) on the main chain.
    # We look at the internal carbons (excluding the acid carbon and the chain terminus).
    branch_found = False
    branch_reason = ""
    max_branch_atoms = 2  # Allowed branch size (1 for methyl, 2 for ethyl)
    
    for atom_idx in main_chain[1:-1]:
        atom = mol.GetAtomWithIdx(atom_idx)
        # Check each carbon neighbor not already on the main chain.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            nbr_idx = nbr.GetIdx()
            if nbr_idx in acyl_chain_set:
                continue
            # Found a candidate branch start. Count connected carbons (restricted to those outside of main chain).
            visited = set()
            stack = [nbr_idx]
            while stack:
                current = stack.pop()
                if current in visited:
                    continue
                visited.add(current)
                current_atom = mol.GetAtomWithIdx(current)
                for nb in current_atom.GetNeighbors():
                    if nb.GetAtomicNum() == 6:
                        nb_idx = nb.GetIdx()
                        if nb_idx in acyl_chain_set or nb_idx in visited:
                            continue
                        stack.append(nb_idx)
            branch_size = len(visited)
            if branch_size <= max_branch_atoms:
                branch_found = True
                branch_reason = f"Found a branch of {branch_size} carbon(s) off the main chain."
                break
        if branch_found:
            break

    if not branch_found:
        return False, "No small alkyl substituents (branch) found on the fatty acyl chain."
    
    return True, ("Contains a carboxylic acid group with a dominant fatty acyl chain "
                  "and a small alkyl branch. " + branch_reason)


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
        ("CCCC\\C(=C/CC)C(O)=O", "2-n-Propyl-2-pentenoic acid"),
        ("CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\C(O)=O", "(2E,6E,10E)-geranylgeranic acid"),
        ("CC(CCCCCCC/C=C/C(=O)O)C", "(E)-11-methyldodec-2-enoic acid"),
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
        ("CC(CO)CCCC(C)CCCC(C)CCCC(C)CC(O)=O", "omega-hydroxyphytanic acid"),
        ("OC(=O)/C=C\\C(C)(C)C", "4,4-dimethyl-2Z-pentenoic acid"),
        ("C([C@H]([C@@H](CCCCCCCCCCCCCCCC[C@@H]1[C@H](CCCCCCCCCCCCCCC[C@@H](O)[C@H](C(O)=O)CCCCCCCCCCCCCCCCCCCCCCCC)C1)OC)C)CCCCCCCCCCCCCCCCC", 
         "(2R)-2-[(1R)-1-hydroxy-16-{(1R,2S)-2-[(17R,18R)-17-methoxy-18-methylhexatriacontyl]cyclopropyl}hexadecyl]hexacosanoic acid"),
        # False positive candidate:
        ("OC(=O)CCC(CCCC)CC", "4-Ethyloctanoic acid"),
        # False negatives:
        ("[H][C@@]12[C@H](CCN1CC=C2COC(=O)[C@](O)([C@H](C)O)C(C)(C)O)OC(=O)C(\\C)=C/C", "heliosupine"),
        ("[O-]C(=O)CCC([N+](C)(C)C)C", "4-aminovaleric acid betaine")
    ]
    
    for smi, name in examples:
        result, reason = is_branched_chain_fatty_acid(smi)
        print(f"SMILES: {smi}\nNAME: {name}\nResult: {result}\nReason: {reason}\n{'-'*70}")