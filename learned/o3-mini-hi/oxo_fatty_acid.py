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
  - It must contain at least one extra carbonyl (ketone or aldehyde) that is NOT part
    of the acid group and that is not located at a terminal carbon of the “fatty acid” chain.
  - In addition, the molecule’s carbon skeleton (the subgraph of carbon atoms) is
    examined. If the longest chain covers fewer than 70% of all carbons, the molecule
    is judged too branched and thus not a typical (linear) fatty acid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import deque

def get_longest_carbon_chain(mol):
    """
    Given an RDKit molecule, constructs a graph over carbon atoms (by index) and
    returns the length of the longest path (in number of atoms) and the endpoints
    (as a set of atom indices) corresponding to one such longest chain.
    (We use a two-pass BFS method which is valid for trees or acyclic graphs.)
    """
    # Get indices of carbon atoms
    carbon_idx = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_idx:
        return 0, set()
    # Build a dict of neighbors for carbons only:
    neighbors = {idx: [] for idx in carbon_idx}
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            neighbors[a1.GetIdx()].append(a2.GetIdx())
            neighbors[a2.GetIdx()].append(a1.GetIdx())
    
    # Pick an arbitrary carbon as starting point.
    start = carbon_idx[0]

    # BFS from start to find the farthest carbon
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
                    queue.append((nb, dist+1, path + [nb]))
        return farthest_node, farthest_dist, farthest_path

    node1, _, _ = bfs(start)
    node2, dist, path = bfs(node1)
    # Return chain length in terms of number of atoms and endpoints (first and last in path)
    endpoints = {path[0], path[-1]}
    return len(path), endpoints

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    
    Heuristics:
      - The molecule must contain only C, H, and O.
      - It should be acyclic.
      - It should have at least 5 carbons.
      - It must have exactly one carboxylic acid group (the terminal -C(=O)[OH]).
      - It must contain at least one extra carbonyl (ketone or aldehyde) that is
        not part of the acid group and that is not at a terminal position in the carbon chain.
      - The carbon skeleton (graph over carbon atoms) should be “linear enough”: i.e.
        the longest carbon chain should cover at least 70% of all carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Must contain only C, H, and O atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6, 8):
            return False, "Molecule contains elements other than C, H, and O."
    
    # Must be acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings and does not appear to be a linear fatty acid."
    
    # Must have at least 5 carbon atoms.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 5:
        return False, f"Molecule contains too few carbons ({carbon_count}) to be a fatty acid."

    # Check carbon skeleton linearity. Compute longest carbon chain.
    longest_chain_len, chain_endpoints = get_longest_carbon_chain(mol)
    if longest_chain_len / carbon_count < 0.70:
        return False, "Molecule appears too branched to be a typical linear fatty acid."

    # Look for carboxylic acid group.
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if len(acid_matches) == 0:
        return False, "Molecule does not contain a carboxylic acid group."
    elif len(acid_matches) > 1:
        return False, "Molecule contains more than one carboxylic acid group, not a typical fatty acid."
    
    # For the single acid group, record the carbonyl atom index (first atom in match).
    acid_carbonyl_idx = acid_matches[0][0]

    # Now search for an extra carbonyl function (ketone or aldehyde) that is not part of the acid
    # and is not at a terminal carbon of the carbon chain.
    extra_oxo_found = False

    # Use SMARTS for ketone:
    ketone_smarts = "[#6][CX3](=O)[#6]"
    ketone_pattern = Chem.MolFromSmarts(ketone_smarts)
    for match in mol.GetSubstructMatches(ketone_pattern):
        # match is a tuple (R, C=O, R'); if the carbonyl carbon (match[1]) is the same as acid carbonyl, skip.
        carbonyl_idx = match[1]
        if carbonyl_idx == acid_carbonyl_idx:
            continue
        # Also require that this extra oxo carbon is not one of the endpoints of the longest chain.
        if carbonyl_idx in chain_endpoints:
            continue
        extra_oxo_found = True
        break

    # If no ketone found, try for aldehyde.
    if not extra_oxo_found:
        aldehyde_smarts = "[#6][CX3H1](=O)"
        aldehyde_pattern = Chem.MolFromSmarts(aldehyde_smarts)
        for match in mol.GetSubstructMatches(aldehyde_pattern):
            carbonyl_idx = match[1]
            if carbonyl_idx == acid_carbonyl_idx:
                continue
            if carbonyl_idx in chain_endpoints:
                continue
            extra_oxo_found = True
            break

    if not extra_oxo_found:
        return False, "Molecule does not appear to contain an additional (internal) aldehydic or ketonic group."

    return True, ("Molecule is an oxo fatty acid: it is acyclic, contains only C, H and O, has a single "
                  "carboxylic acid group and an extra (non‐terminal) oxo function.")

# (Optional testing block)
if __name__ == '__main__':
    test_smiles = [
        "CC(C)C(=O)CC(O)=O",  # 4-methyl-3-oxopentanoic acid (true)
        "O=C(CCCC)CCC(O)=O",   # 4-keto-n-caprylic acid (false by our extra criteria; extra oxo at chain end)
        "CCCCCCCCCCCCC(=O)CC(O)=O",  # 3-oxopalmitic acid (true)
        "OC(=O)CCCCCCCCCCCCCCCCCCCCC=O",  # 22-oxodocosanoic acid (true)
        "O=C(CCC(O)=O)CC",     # 4-keto-n-caproic acid (false extra oxo terminal)
    ]
    for smi in test_smiles:
        result, reason = is_oxo_fatty_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")