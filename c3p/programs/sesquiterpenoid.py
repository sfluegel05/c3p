"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
"""
Classifies: CHEBI:35189 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    Sesquiterpenoids are terpenoids derived from sesquiterpenes (C15 skeleton),
    possibly modified by rearrangement or removal of some carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Sesquiterpenoids should have around 15 carbons (allowing variation due to modifications)
    if c_count < 11 or c_count > 25:
        return False, f"Carbon count ({c_count}) outside typical range for sesquiterpenoids (11-25)"

    # Calculate rings
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    
    # Most sesquiterpenoids are cyclic, but some can be linear
    if ring_count > 5:
        return False, f"Too many rings ({ring_count}) for a typical sesquiterpenoid"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 180 or mol_wt > 600:
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical range for sesquiterpenoids"

    # Count double bonds using SMARTS pattern
    double_bond_pattern = Chem.MolFromSmarts('C=C')
    if double_bond_pattern:
        double_bond_count = len(mol.GetSubstructMatches(double_bond_pattern))
    else:
        double_bond_count = 0

    # Look for typical functional groups found in sesquiterpenoids
    functional_groups = {
        'alcohol': '[OH1]',
        'ketone': '[CX3](=[OX1])[CX4]',
        'carboxylic_acid': '[CX3](=O)[OX2H1]',
        'ester': '[#6][CX3](=O)[OX2H0][#6]',
        'epoxide': '[C]1O[C]1',
        'ether': '[OX2]([CX4])[CX4]'
    }
    
    found_groups = []
    for group_name, smarts in functional_groups.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            found_groups.append(group_name)

    # Look for branching methyl groups - common in terpenoids
    methyl_pattern = Chem.MolFromSmarts('[CH3]')
    if methyl_pattern:
        methyl_count = len(mol.GetSubstructMatches(methyl_pattern))
    else:
        methyl_count = 0

    # Common sesquiterpenoid ring patterns
    ring_patterns = [
        ('[C]1[C][C][C][C][C]1', 'cyclohexane'),  # 6-membered ring
        ('[C]1[C][C][C][C]1', 'cyclopentane'),    # 5-membered ring
        ('[C]1[C][C][C]1', 'cyclobutane'),        # 4-membered ring
    ]
    
    found_rings = []
    for smarts, ring_name in ring_patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            found_rings.append(ring_name)

    # Scoring system
    score = 0
    reasons = []
    
    # Base score from carbon count
    if 13 <= c_count <= 17:
        score += 2
        reasons.append(f"Typical sesquiterpenoid carbon count ({c_count})")
    elif 11 <= c_count <= 25:
        score += 1
        reasons.append(f"Acceptable carbon count ({c_count})")

    # Score from rings
    if 1 <= ring_count <= 4:
        score += 2
        reasons.append(f"Typical ring count ({ring_count})")
    elif ring_count == 0 and double_bond_count >= 2:
        score += 1
        reasons.append("Linear with multiple double bonds")

    # Score from functional groups
    if found_groups:
        score += len(found_groups)
        reasons.append(f"Contains {', '.join(found_groups)}")

    # Score from methyl groups
    if 3 <= methyl_count <= 6:
        score += 2
        reasons.append(f"Contains {methyl_count} methyl groups")
    elif methyl_count > 0:
        score += 1
        reasons.append(f"Contains {methyl_count} methyl groups")

    # Final decision
    if score >= 4:
        return True, "; ".join(reasons)
    else:
        return False, "Insufficient sesquiterpenoid characteristics"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26658',
                          'name': 'sesquiterpenoid',
                          'definition': 'Any terpenoid derived from a '
                                        'sesquiterpene. The term includes '
                                        'compounds in which the C15 skeleton '
                                        'of the parent sesquiterpene has been '
                                        'rearranged or modified by the removal '
                                        'of one or more skeletal atoms '
                                        '(generally methyl groups).',
                          'parents': ['CHEBI:26873'],
                          'xrefs': ['LIPID_MAPS_class:LMPR0103'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               "Error: module 'rdkit.Chem.rdMolDescriptors' has no attribute "
               "'CalcNumDoubleBonds'\n"
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: NONE\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'O=C(O)[C@]12C(=C[C@@H]3[C@]1(C[C@H]4[C@H]([C@@]2(C=O)C3)CC[C@H]4C)CO[C@]56O[C@@H]7O[C@H]([C@@H]5O)[C@@H]([C@@H]7O6)OC(=O)C)C(C)C',
                                      'name': 'BE-31405',
                                      'reason': 'Carbon count (28) outside '
                                                'typical range for '
                                                'sesquiterpenoids (11-25)'},
                                  {   'smiles': 'O=C1C2=C([C@@]3(CC[C@H]4[C@]3(C1)CO[C@@]4(CCC(=O)CC)C)C)[C@@H](OC(=O)C)C[C@H]([C@@]2(CCC(=O)O)C)C(=C)C',
                                      'name': 'Ganoboninone H',
                                      'reason': 'Carbon count (31) outside '
                                                'typical range for '
                                                'sesquiterpenoids (11-25)'},
                                  {   'smiles': 'O=C1C(=O)C(=C(O)C(=C1/C=C/C2=CC=CC=C2)C3=C(O)C(=C(O)C=C3/C=C/C4=CC=CC=C4)C(C)C)C(C)C',
                                      'name': 'Duotap-520',
                                      'reason': 'Carbon count (34) outside '
                                                'typical range for '
                                                'sesquiterpenoids (11-25)'},
                                  {   'smiles': '[H][C@]12C[C@H](OC(=O)c3ccccc3)[C@]3(C)[C@@H](OC(C)=O)[C@H](C[C@@H](C)[C@@]3(OC1(C)C)[C@@H]2OC(C)=O)OC(C)=O',
                                      'name': 'Triptogelin C-1',
                                      'reason': 'Carbon count (28) outside '
                                                'typical range for '
                                                'sesquiterpenoids (11-25)'},
                                  {   'smiles': 'O=C1[C@H]([C@]2([C@@]([C@](CC2)([C@H](C)/C=C/CC(C)C)[H])(CC1)C)[H])CCC3=C(C=CC(O)=C3)C',
                                      'name': 'calicoferol A',
                                      'reason': 'Carbon count (27) outside '
                                                'typical range for '
                                                'sesquiterpenoids (11-25)'},
                                  {   'smiles': 'ClC1=C2N([C@H]3N[C@H](C(=O)O[C@@H]4[C@]5(O)[C@H](C=C([C@H]4OC(=O)C)C)[C@H](C(=C)C)CC[C@H]5C)C[C@]3(C2=CC=C1)O)C',
                                      'name': '7-chlorofischerindoline',
                                      'reason': 'Carbon count (29) outside '
                                                'typical range for '
                                                'sesquiterpenoids (11-25)'},
                                  {   'smiles': 'O=C1C(=C2[C@@](CC2)(C)[C@@H]3[C@@]1(O)CC(C)(C)C3)COC(=O)CCCCCCC/C=C/CCCCCCCC',
                                      'name': 'Atlanticone A',
                                      'reason': 'Carbon count (33) outside '
                                                'typical range for '
                                                'sesquiterpenoids (11-25)'},
                                  {   'smiles': 'O=C1NCCCC2NC(=O)C(C2=O)=C(O)C=CC3C(CC=C1)C4C5C(C(C(C(OC)C)C(C5)C)C6C4O6)C3',
                                      'name': 'Ripromycin',
                                      'reason': 'Carbon count (30) outside '
                                                'typical range for '
                                                'sesquiterpenoids (11-25)'},
                                  {   'smiles': 'O=C1NC(OCCCC)(C)C=C1C(=O)[C@@H]2[C@@H](C(=C[C@H]3[C@@H]2[C@H](C[C@H](C3)C)C)C)/C(=C/C)/C',
                                      'name': 'Ascosalipyrrolidinone A',
                                      'reason': 'Carbon count (27) outside '
                                                'typical range for '
                                                'sesquiterpenoids (11-25)'},
                                  {   'smiles': 'O=C1OCC23C4(C5(OC5)C(OC2C=C(C)CC3)C(C4OC(=O)C=CC=CC6OC(CC(=C1)C)OC6C)O)C',
                                      'name': 'Roridin M',
                                      'reason': 'Carbon count (29) outside '
                                                'typical range for '
                                                'sesquiterpenoids (11-25)'}],
    'attempt': 1,
    'reasoning': 'Let me analyze the previous attempt and suggest '
                 'improvements:\n'
                 '\n'
                 '1. The main error was using `CalcNumDoubleBonds` which '
                 "doesn't exist in rdMolDescriptors. We can instead use a "
                 'SMARTS pattern to count double bonds.\n'
                 '\n'
                 '2. Looking at the example structures provided:\n'
                 '- They all have a C15 base skeleton (some modified)\n'
                 '- Many contain oxygen-containing functional groups '
                 '(alcohols, ketones, esters)\n'
                 '- Most are cyclic with 1-4 rings\n'
                 '- They have varying degrees of unsaturation\n'
                 '- Some have modifications that add or remove carbons from '
                 'the base C15 skeleton\n'
                 '\n'
                 '3. Improvements needed:\n'
                 '- Fix the double bond counting\n'
                 '- Adjust carbon count range (some examples have more carbons '
                 'due to modifications)\n'
                 '- Add more specific structural patterns for sesquiterpenoid '
                 'skeletons\n'
                 '- Consider molecular weight range adjustment\n'
                 '- Add checks for specific ring systems common in '
                 'sesquiterpenoids\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 558,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.8480243161094225,
    'f1': 0.9177631578947368,
    'accuracy': 0.8480243161094225,
    'negative_predictive_value': 0.0}