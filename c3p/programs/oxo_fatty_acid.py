"""
Classifies: CHEBI:59644 oxo fatty acid
"""
"""
Classifies: Oxo fatty acid
Definition: Any fatty acid containing at least one aldehydic or ketonic group,
            in addition to the carboxylic acid group.
Heuristic criteria used here:
  - The molecule must contain only C, H, and O.
  - It must be acyclic.
  - It must have at least a minimal number of carbons (>=5).
  - It must have exactly one carboxylic acid group (SMARTS "[CX3](=O)[OX2H]").
  - It must contain at least one extra carbonyl (ketone or aldehyde) that is not
    the carboxylic acid.
  - If the extra carbonyl function is located at a terminal position of the longest
    carbon chain, then the chain is required to be long (≥10 carbons) to avoid false positives.
  - The carbon skeleton (graph over carbon atoms) should be “linear enough”:
    the longest carbon chain must cover at least 70% of all carbons.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import deque

def get_longest_carbon_chain(mol):
    """
    Given an RDKit molecule, builds a graph over carbon atoms and returns the length (number of atoms)
    of the longest chain and the endpoints (as a set of atom indices) of one such chain.
    """
    # Get all carbon atom indices
    carbon_idx = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_idx:
        return 0, set()
    # Build neighbor dict for carbons only
    neighbors = {idx: [] for idx in carbon_idx}
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            neighbors[a1.GetIdx()].append(a2.GetIdx())
            neighbors[a2.GetIdx()].append(a1.GetIdx())
    
    # Helper BFS function starting from given index
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

    start = carbon_idx[0]
    node1, _, _ = bfs(start)
    node2, dist, path = bfs(node1)
    endpoints = {path[0], path[-1]}
    return len(path), endpoints

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    
    Heuristics:
      - The molecule must consist only of C, H, and O atoms.
      - It must be acyclic.
      - It must have at least 5 carbon atoms.
      - Its carbon skeleton must be linear enough: the longest carbon chain must
        cover at least 70% of all carbon atoms.
      - There must be exactly one carboxylic acid group (SMARTS "[CX3](=O)[OX2H]").
      - There must be at least one extra carbonyl (ketone or aldehyde) that is not the acid function.
      - If an extra carbonyl falls at a terminal position of the carbon chain (i.e. one of the chain endpoints)
        and is not the acid carbonyl, we only accept it if the longest chain is sufficiently long (>=10 carbons).
    
    Returns:
        (bool, str): Tuple containing classification result and reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # 1. Check that the molecule contains only C, H, and O.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6, 8):
            return False, "Molecule contains elements other than C, H, and O."
    
    # 2. Must be acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings and does not appear to be a linear fatty acid."
    
    # 3. Must have at least 5 carbons.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 5:
        return False, f"Molecule contains too few carbons ({carbon_count}) to be a fatty acid."
    
    # 4. Check overall linearity using the longest carbon chain.
    longest_chain_len, chain_endpoints = get_longest_carbon_chain(mol)
    if longest_chain_len / carbon_count < 0.70:
        return False, "Molecule appears too branched to be a typical linear fatty acid."

    # 5. Search for exactly one carboxylic acid group.
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if len(acid_matches) == 0:
        return False, "Molecule does not contain a carboxylic acid group."
    elif len(acid_matches) > 1:
        return False, "Molecule contains more than one carboxylic acid group, not a typical fatty acid."
    
    # Record the carbonyl atom of the acid group.
    acid_carbonyl_idx = acid_matches[0][0]

    # 6. Search for an extra carbonyl function (ketone or aldehyde) that is not the acid group.
    extra_oxo_found = False
    extra_reason = ""
    accepted_candidate = None  # store candidate candidate if found

    # Define patterns for ketone and aldehyde.
    ketone_smarts = "[#6][CX3](=O)[#6]"
    ketone_pattern = Chem.MolFromSmarts(ketone_smarts)
    aldehyde_smarts = "[#6][CX3H1](=O)"
    aldehyde_pattern = Chem.MolFromSmarts(aldehyde_smarts)
    
    # Helper routine to check candidate carbonyl:
    def candidate_is_valid(carbonyl_idx):
        # Skip if carbonyl is the acid carbonyl.
        if carbonyl_idx == acid_carbonyl_idx:
            return False
        # If the candidate is at a chain endpoint, then require a sufficiently long chain.
        if carbonyl_idx in chain_endpoints:
            if longest_chain_len < 10:
                return False
        return True

    # Check ketone candidates.
    for match in mol.GetSubstructMatches(ketone_pattern):
        # match: (R, C(=O), R'). Examine the carbonyl atom (middle).
        carbonyl_idx = match[1]
        if candidate_is_valid(carbonyl_idx):
            extra_oxo_found = True
            accepted_candidate = carbonyl_idx
            break

    # If no valid ketone found, try for aldehyde.
    if not extra_oxo_found:
        for match in mol.GetSubstructMatches(aldehyde_pattern):
            # match: (R, C(=O)H). The carbonyl atom is in position 0 of the =O?
            # In this SMARTS, the carbonyl carbon is the first atom in the pattern.
            carbonyl_idx = match[0]
            if candidate_is_valid(carbonyl_idx):
                extra_oxo_found = True
                accepted_candidate = carbonyl_idx
                break

    if not extra_oxo_found:
        return False, "Molecule does not appear to contain an additional (non–acid) aldehydic or ketonic group."
    
    return True, ("Molecule is an oxo fatty acid: it is acyclic, contains only C, H and O, has a single "
                  "carboxylic acid group and an extra (aldehydic or ketonic) function.")

# (Optional testing block)
if __name__ == '__main__':
    test_smiles = [
        # True positives:
        "CC(C)C(=O)CC(O)=O",  # 4-methyl-3-oxopentanoic acid
        "CCCC(O)C(=O)CCC(=O)C\\C=C/CCCCCC(O)=O",  # (7Z)-14-hydroxy-10,13-dioxoheptadec-7-enoic acid
        "CCCCCCCCCCCCC(=O)CC(O)=O",  # 3-oxopalmitic acid
        "C(O)(CCC(CCCC\\C=C/C=C\\C=C\\C=C\\CC)=O)=O",  # (9Z,11Z,13E,15E)-4-oxooctadecatetraenoic acid
        "CCC(O)CC(=O)C(O)=O",  # 4-hydroxy-2-oxohexanoic acid
        "CCCCCC\\C=C/CC(=O)\\C=C\\C=C/C\\C=C/CCCC(O)=O",  # 12-oxo-ETE
        "CCCC(=O)C(O)=O",  # 2-oxopentanoic acid
        "CCC(C)(O)C(=O)C(O)=O",  # 3-hydroxy-3-methyl-2-oxopentanoic acid
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
        "C(C(O)=O)CCCC/C=C\\C/C=C\\CC(/C=C/C=C\\C/C=C\\CC)=O",  # (7Z,10Z,14E,16Z,19Z)-13-oxodocosapentaenoic acid
        "C(CCCCCCC(CC(=O)O)=O)CCCCCCCCCCCCCCCC",  # 3-oxohexacosanoic acid
        "C(C(O)=O)C/C=C\\C/C=C\\C/C=C\\C\\C=C/C=C/C(C/C=C\\CC)=O",  # (4Z,7Z,10Z,13Z,15E,19Z)-17-oxodocosahexaenoic acid
        # False positives examples (expected False)
        "C(CCCCCC(=O)O)C/C=C/CC(/C=C(/C=C/CCO)\\O)=O",  # (9E,13Z,15E)-14,18-dihydroxy-12-oxooctadeca-9,13,15-trienoic acid
        "O=C(CCCC)CCC(O)=O",  # 4-keto-n-caprylic acid
        "O=C(CCC(O)=O)CC",  # 4-keto-n-caproic acid
    ]
    for smi in test_smiles:
        result, reason = is_oxo_fatty_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")