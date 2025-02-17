"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
"""
Classifies: branched-chain fatty acid (BCFA)
Definition:
  A branched-chain fatty acid is a fatty acid – that is, a molecule having
  a carboxylic acid functional group plus a long alkyl (fatty acyl) chain –
  in which the fatty acyl chain carries one or more small extra alkyl substituents.
  Typically the substituent is a one‐carbon (methyl) group, though in some cases
  a two‐carbon (ethyl) group is observed.
  
This implementation works as follows:
  1. It first searches for a carboxylic acid. In addition to standard patterns,
     we include an extra SMARTS to catch representations where the –OH is drawn 
     as connected via an oxygen to the acid carbon.
  2. Once an acid group is located (using one of the SMARTS), we select the acid carbon.
  3. We then construct a carbon-only graph of the molecule and – starting from any carbon
     attached to the acid carbon – we enumerate all simple (acyclic) paths (via DFS)
     that reach a “terminal” carbon (one with only one carbon neighbour in this graph).
     Of these, we select the longest path and treat that as the dominant fatty acyl chain.
  4. Lastly, for each internal carbon along that chain (i.e. excluding the first and last atoms)
     we look for substituents (neighbors not on the main chain) and compute their connected 
     carbon size (using DFS that stays off the main chain). If we find at least one branch 
     of size 1 or 2 (i.e. methyl or ethyl), we classify it as a branched-chain fatty acid.

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
        bool: True if molecule qualifies as a BCFA, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
        
    # Step 1. Look for a carboxylic acid group. We use several SMARTS patterns 
    # to catch both standard and alternate representations.
    acid_smarts_list = [
        "[CX3](=O)[OX2H]",    # protonated acid: C(=O)OH
        "[CX3](=O)[OX2-]",    # deprotonated acid: C(=O)[O-]
        "[OX2H][CX3](=O)[#6]"  # sometimes the -OH is drawn attached to O which in turn attaches to C(=O)C...
    ]
    acid_matches = []
    acid_pattern_used = None
    for smt in acid_smarts_list:
        patt = Chem.MolFromSmarts(smt)
        matches = mol.GetSubstructMatches(patt)
        if matches:
            acid_matches = matches
            acid_pattern_used = smt
            break
    if not acid_matches:
        return False, "No carboxylic acid group found; not a fatty acid."
    
    # We assume that in our acid SMARTS the acid carbon is atom position 0 
    # (or in the alternative pattern we choose the [CX3](=O) part).
    acid_idx = acid_matches[0][0]
    
    # Step 2. Build a carbon-only graph.
    all_carbons = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(all_carbons) < 3:
        return False, "Too few carbon atoms to be a fatty acid."
    
    # Create a dictionary mapping each carbon index to its carbon neighbors.
    carbon_graph = { idx: [] for idx in all_carbons }
    for idx in all_carbons:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                carbon_graph[idx].append(nbr.GetIdx())
    
    # Step 3. Identify fatty acyl chain.
    # The acyl chain should start from a carbon neighbor of the acid carbon.
    acid_neighbors = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(acid_idx).GetNeighbors() if nbr.GetAtomicNum() == 6]
    if not acid_neighbors:
        return False, "Acid carbon not attached to any carbon, so not a fatty acid."
    
    # Define "terminal" carbons in the carbon_graph as those with only 1 carbon neighbor
    # (excluding the acid carbon if it is included).
    terminal_carbons = set(idx for idx in all_carbons if len(carbon_graph[idx]) == 1 and idx != acid_idx)
    
    # DFS function to get all simple paths from a start carbon to any terminal carbon.
    def dfs_paths(current, path, visited):
        paths = []
        # if we are at a terminal (and not the very beginning)
        if current in terminal_carbons and len(path) > 1:
            paths.append(path)
            # we do not return immediately because other longer paths might exist.
        for nbr in carbon_graph.get(current, []):
            if nbr not in visited:
                new_visited = visited | {nbr}
                paths.extend(dfs_paths(nbr, path + [nbr], new_visited))
        return paths

    candidate_chains = []
    # For each acid neighbor, compute the simple paths that reach a terminal carbon.
    for start in acid_neighbors:
        paths = dfs_paths(start, [start], {acid_idx, start})
        if paths:
            # choose the longest path from this start
            longest = max(paths, key=lambda p: len(p))
            candidate_chains.append(longest)
    if not candidate_chains:
        return False, "Could not extract a fatty acyl chain from the acid carbon."
    
    # Choose the candidate chain with the maximum length.
    main_chain = max(candidate_chains, key=lambda p: len(p))
    # For checking branches later, use a set for carbons in the acyl chain.
    acyl_chain_set = set(main_chain)
    
    # Step 3a. Check that the main chain is long enough (expect at least 3 carbons)
    if len(main_chain) < 3:
        return False, "The fatty acyl chain appears too short."
    
    # Step 3b. Check that the main chain is dominant if the molecule is large.
    ratio = len(acyl_chain_set) / len(all_carbons)
    if len(all_carbons) > 8 and ratio < 0.5:
        return False, "The fatty acyl chain does not dominate the carbon skeleton."
        
    # Step 4. Look for a small alkyl branch on the main chain.
    # We examine the internal carbons (excluding the first and the terminal carbon).
    branch_found = False
    branch_reason = ""
    max_branch_atoms = 2  # allowed branch: methyl (1) or ethyl (2)

    # For each internal atom in the acyl chain, check for extra carbon neighbors not in the chain.
    for atom_idx in main_chain[1:-1]:
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            nbr_idx = nbr.GetIdx()
            if nbr_idx in acyl_chain_set:
                continue  # this neighbor is part of the main chain
            # Found a candidate branch start.
            # Using DFS (restricted to the branch not merging back with the main chain)
            visited = set()
            stack = [nbr_idx]
            while stack:
                current = stack.pop()
                if current in visited:
                    continue
                visited.add(current)
                for nb in carbon_graph.get(current, []):
                    if nb in acyl_chain_set or nb in visited:
                        continue
                    stack.append(nb)
            branch_size = len(visited)
            if branch_size <= max_branch_atoms:
                branch_found = True
                branch_reason = f"Found a branch of {branch_size} carbon(s) off the main chain."
                break
        if branch_found:
            break

    if not branch_found:
        return False, "No small alkyl substituents (branch) found on the fatty acyl chain."
    
    # All checks passed.
    explanation = ("Contains a carboxylic acid group with a dominant fatty acyl chain and a small alkyl branch. "
                   + branch_reason)
    return True, explanation


# Example usage (for quick testing)
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
        # Some false negatives:
        ("[H][C@@]12[C@H](CCN1CC=C2COC(=O)[C@](O)([C@H](C)O)C(C)(C)O)OC(=O)C(\\C)=C/C", "heliosupine"),
        ("[O-]C(=O)CCC([N+](C)(C)C)C", "4-aminovaleric acid betaine")
    ]
    
    for smi, name in examples:
        result, reason = is_branched_chain_fatty_acid(smi)
        print(f"SMILES: {smi}\nNAME: {name}\nResult: {result}\nReason: {reason}\n{'-'*70}")