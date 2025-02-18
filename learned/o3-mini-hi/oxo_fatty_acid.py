"""
Classifies: CHEBI:59644 oxo fatty acid
"""
"""
Classifies: Oxo fatty acid
Definition: Any fatty acid containing at least one aldehydic or ketonic group 
            in addition to the carboxylic acid group.
Heuristic criteria used here (improved):
  - Molecule must consist only of C, H, and O.
  - Molecule must be acyclic.
  - There must be at least 5 carbon atoms.
  - The overall carbon skeleton should be “linear enough”: 
      the longest carbon chain must span at least 70% of all carbons.
  - The molecule must contain exactly one carboxylic acid group (SMARTS "[CX3](=O)[OX2H]"),
      and its acid carbon must be an endpoint of the longest chain.
  - There must be at least one additional carbonyl (ketonic or aldehydic) function
      (SMARTS "[#6][CX3](=O)[#6]" or "[#6][CX3H1](=O)") that is not the acid carbon.
      Furthermore, the extra oxo must not lie too “close” (topologically) to the acid (min. path length 3),
      and if it is at a chain endpoint then the longest chain must be at least 10 carbons long.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import deque

def get_longest_carbon_chain(mol):
    """
    Build a graph over carbon atoms and calculate the longest chain (path)
    in terms of number of carbon atoms and return the length and endpoints.
    """
    # Gather indices for all carbons.
    carbon_idx = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_idx:
        return 0, set()
    # Build neighbor dictionary for carbons only.
    neighbors = {idx: [] for idx in carbon_idx}
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            neighbors[a1.GetIdx()].append(a2.GetIdx())
            neighbors[a2.GetIdx()].append(a1.GetIdx())
    
    # For each carbon, do a simple BFS to estimate the farthest carbon.
    def bfs(start_idx):
        visited = {start_idx}
        queue = deque([(start_idx, 0, [start_idx])])
        farthest_node = start_idx
        farthest_dist = 0
        farthest_path = [start_idx]
        while queue:
            node, dist, path = queue.popleft()
            if dist > farthest_dist:
                farthest_dist = dist
                farthest_node = node
                farthest_path = path
            for nb in neighbors[node]:
                if nb not in visited:
                    visited.add(nb)
                    queue.append((nb, dist + 1, path + [nb]))
        return farthest_node, farthest_dist, farthest_path
    
    # Start from an arbitrary carbon.
    start = carbon_idx[0]
    node1, _, _ = bfs(start)
    node2, dist, path = bfs(node1)
    endpoints = {path[0], path[-1]}
    return len(path), endpoints

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    
    Heuristics (improved):
      - Molecule must contain only C, H, and O.
      - Molecule must be acyclic.
      - Must have at least 5 carbon atoms.
      - Its carbon skeleton must be “linear enough”:
           longest carbon chain must cover at least 70% of all carbons.
      - There must be exactly one carboxylic acid group (SMARTS "[CX3](=O)[OX2H]").
        Furthermore, the acid group’s carbonyl must be at an endpoint of the longest chain.
      - There must be at least one extra carbonyl (ketone or aldehyde) that is not part of the acid.
         For any candidate extra oxo, we require that:
           * Its topological distance (number of bonds) from the acid carbon is >= 3.
           * If the candidate lies at a chain endpoint, then the longest chain must have at least 10 carbons.
    
    Returns:
        (bool, str): Tuple with classification result and reasoning.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string."
    
    # (1) Check that atoms are only C, H, and O.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6, 8):
            return False, "Molecule contains elements other than C, H, and O."
    
    # (2) Require acyclic structure.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings and does not appear to be a linear fatty acid."
    
    # (3) Must have a sufficient number of carbons.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 5:
        return False, f"Molecule contains too few carbons ({carbon_count}) to be a fatty acid."
    
    # (4) Check overall linearity: longest carbon chain covers at least 70% of all carbons.
    longest_chain_len, chain_endpoints = get_longest_carbon_chain(mol)
    if longest_chain_len / carbon_count < 0.70:
        return False, "Molecule appears too branched to be a typical linear fatty acid."
    
    # (5) Identify carboxylic acid group (exactly one).
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if len(acid_matches) == 0:
        return False, "Molecule does not contain a carboxylic acid group."
    if len(acid_matches) > 1:
        return False, "Molecule contains more than one carboxylic acid group."
    
    # We assume the first matching acid group.
    # The acid carbon is taken as the carbonyl carbon (first atom in the match).
    acid_carbonyl_idx = acid_matches[0][0]
    if acid_carbonyl_idx not in chain_endpoints:
        return False, "The carboxylic acid group is not located at the terminus of the main carbon chain."
    
    # (6) Look for an additional carbonyl (ketone or aldehyde) that is not the acid.
    # Define SMARTS for ketone and aldehyde.
    ketone_smarts = "[#6][CX3](=O)[#6]"
    aldehyde_smarts = "[#6][CX3H1](=O)"
    ketone_pattern = Chem.MolFromSmarts(ketone_smarts)
    aldehyde_pattern = Chem.MolFromSmarts(aldehyde_smarts)
    
    # Helper: require extra carbonyl is not the acid carbon,
    # its shortest path from the acid carbon is at least 3 bonds,
    # and if it lies at an endpoint of the longest chain, then the chain must be long (>=10).
    def candidate_is_valid(carbonyl_idx):
        if carbonyl_idx == acid_carbonyl_idx:
            return False
        # Compute topological distance (number of bonds) between acid carbon and candidate.
        path = Chem.GetShortestPath(mol, acid_carbonyl_idx, carbonyl_idx)
        if len(path) < 3:  # path length counts atoms along the shortest route (so <3 means less than 2 bonds apart)
            return False
        if carbonyl_idx in chain_endpoints and longest_chain_len < 10:
            return False
        return True

    extra_oxo_found = False
    # Try ketone candidates.
    for match in mol.GetSubstructMatches(ketone_pattern):
        # In the ketone pattern the carbonyl is the second atom.
        candidate_idx = match[1]
        if candidate_is_valid(candidate_idx):
            extra_oxo_found = True
            break

    # If no valid ketone, try aldehyde candidates.
    if not extra_oxo_found:
        for match in mol.GetSubstructMatches(aldehyde_pattern):
            candidate_idx = match[0]  # in the aldehyde pattern the carbonyl is the first atom.
            if candidate_is_valid(candidate_idx):
                extra_oxo_found = True
                break

    if not extra_oxo_found:
        return False, "Molecule does not appear to contain an additional (non–acid) aldehydic or ketonic group."
    
    return True, ("Molecule is an oxo fatty acid: it is acyclic, contains only C, H and O, "
                  "has a single carboxylic acid group placed at the end of its carbon skeleton, "
                  "and contains an extra (aldehydic or ketonic) function at a suitable location.")

# (Optional testing block)
if __name__ == '__main__':
    test_smiles = [
        # True positives:
        "CC(C)C(=O)CC(O)=O",  # 4-methyl-3-oxopentanoic acid
        "CCCC(O)C(=O)CCC(=O)C\\C=C/CCCCCC(O)=O",  # (7Z)-14-hydroxy-10,13-dioxoheptadec-7-enoic acid
        "CCCCCCCCCCCCC(=O)CC(O)=O",  # 3-oxopalmitic acid
        "OC(=O)CCCCCCCCCCCCCCCCCCCCC=O",  # 22-oxodocosanoic acid
        "C(O)(CCC(CCCC\\C=C/C=C\\C=C\\C=C\\CC)=O)=O",  # (9Z,11Z,13E,15E)-4-oxooctadecatetraenoic acid
        "CCC(O)CC(=O)C(O)=O",       # 4-hydroxy-2-oxohexanoic acid
        "CCCCCC\\C=C/CC(=O)\\C=C\\C=C/C\\C=C/CCCC(O)=O",  # 12-oxo-ETE
        "CCCC(=O)C(O)=O",          # 2-oxopentanoic acid
        "CCC(C)(O)C(=O)C(O)=O",     # 3-hydroxy-3-methyl-2-oxopentanoic acid
        "CCCCCC\\C=C/C\\C=C/C\\C=C/C=C/C(=O)CCCC(O)=O",  # 5-oxo-ETE
        "C(C(/C=C/C=C/C=C\\[C@H](CCCC(O)=O)O)=O)/C=C\\CCCCC=O",  # 12,20-dioxoleukotriene B4
        "C(C(O)=O)CCCCC(/C=C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O",  # (8E,10Z,13Z,16Z,19Z)-7-oxodocosapentaenoic acid
        "C(CCC/C=C\\C/C=C\\CC(/C=C/C=C\\CCCCC)=O)(=O)O",  # 11-oxo-ETE
        "C(CCCCCCCC)CCCCCCC(CC(O)=O)=O",  # 3-oxooctadecanoic acid
        "CC(C)=CCCC(=O)CC(O)=O",  # 7-methyl-3-oxooct-6-enoic acid
        "C(CCCCCCCC(C(=O)O)=O)CC",  # 2-oxododecanoic acid
        "C(C(/C=C/C=C/C=C/[C@H](CCCC(O)=O)O)=O)/C=C\\CCCCC",  # 12-oxo-6-trans-leukotriene B4
        "CC(C)CCCC(C)CCCC(C)CCC(=O)C(C)C(O)=O",  # 3-oxopristanic acid
        "CCCCCC(=O)\\C=C\\C=C/C=C/C=C/[C@@H](O)[C@@H](O)CCCC(O)=O",  # 15-oxolipoxin A4
        "CC[C@@](C)(O)C(=O)C(O)=O",  # (R)-3-hydroxy-3-methyl-2-oxopentanoic acid
        "O=C(CC(=O)O)CC/C=C/[C@@H](C/C=C\\CCCCC)O",  # (6E,8R,10Z)-8-hydroxy-3-oxohexadeca-6,10-dienoic acid
        "C(CCCC(CCCCCCO)=O)CCCCC(O)=O",  # 16-hydroxy-10-oxohexadecanoic acid
        "CCCCCCCCCCCCCCC(=O)C(O)=O",  # 2-oxohexadecanoic acid
        # (other examples omitted for brevity)
        # False positives examples (expected False):
        "C(CCCCCC(=O)O)C/C=C/CC(/C=C(/C=C/CCO)\\O)=O",  # wrongly classified in previous version
        "O=C(CCCC)CCC(O)=O",  # 4-keto-n-caprylic acid
        "O=C(CCC(O)=O)CC",  # 4-keto-n-caproic acid
    ]
    for smi in test_smiles:
        result, reason = is_oxo_fatty_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")