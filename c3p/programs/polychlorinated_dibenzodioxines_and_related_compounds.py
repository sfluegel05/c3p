"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Classifies a molecule as a polychlorinated dibenzodioxin or structurally related compound based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule belongs to the class, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sufficient number of chlorine or bromine atoms
    halogen_count = sum(atom.GetSymbol() in ['Cl', 'Br'] for atom in mol.GetAtoms())
    if halogen_count < 3:  # Lowered to capture more true positives
        return False, "Insufficient halogen atoms for a likely related compound."

    # SMARTS pattern for polychlorinated biphenyls (PCBs)
    biphenyl_pattern = Chem.MolFromSmarts('c1cc(c(cc1)-c2ccccc2)[Cl,Br]')
    
    # SMARTS patterns for polychlorinated dibenzodioxins (PCDDs)
    dioxin_patterns = [
        Chem.MolFromSmarts('Oc1cc2Oc3ccccc3Oc2cc1'),  # Allow variability in Cl positions
        Chem.MolFromSmarts('c1cc2Oc3ccc(Cl)cc3cc2cc1')
    ]
    
    # SMARTS patterns for polychlorinated dibenzofurans (PCDFs)
    furan_patterns = [
        Chem.MolFromSmarts('c1cc2oc3cc(Cl)ccc3c2cc1'),
        Chem.MolFromSmarts('c1cc2oc3ccc(Cl)c2cc3cc1')
    ]
    
    # Check for presence of these patterns
    matches_dioxin_furan = any(mol.HasSubstructMatch(pattern) for pattern in dioxin_patterns + furan_patterns)
    matches_biphenyl = mol.HasSubstructMatch(biphenyl_pattern)

    # Exclude common false patterns
    common_false_patterns = [
        Chem.MolFromSmarts('[OH]'),  # Hydroxy groups are not typical in PCBs, PCDD/Fs
    ]
    false_positive_match = any(mol.HasSubstructMatch(pattern) for pattern in common_false_patterns)

    if (matches_dioxin_furan or matches_biphenyl) and not false_positive_match:
        return True, "Matches polychlorinated dibenzodioxin or related compound pattern."

    return False, "Does not match polychlorinated dibenzodioxin or related compound pattern."


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:134045',
                          'name': 'polychlorinated dibenzodioxines and related '
                                  'compounds',
                          'definition': 'Organochlorine compounds that are '
                                        'polychlorinated dibenzodioxines and '
                                        'structurally related entities that '
                                        'are persistant organic pollutants. '
                                        'These include polychlorinated '
                                        'dibenzofurans as well as '
                                        'polychlorinated and polybrominated '
                                        'biphenyls  They vary widely in their '
                                        'toxicity, but their toxic mode of '
                                        'action is through the aryl '
                                        'hydrocarbon receptor.',
                          'parents': ['CHEBI:17792'],
                          'xrefs': [   'Wikipedia:Dioxins_and_dioxin-like_compounds'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.55 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'Clc1cc2Oc3c(Cl)c(Cl)c(Cl)c(Cl)c3Oc2c(Cl)c1Cl NAME: '
               '1,2,3,4,6,7,8-Heptachlorodibenzodioxin REASON: CORRECT Matches '
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1c(Cl)c(Cl)c2Oc3ccccc3Oc2c1Cl NAME: '
               '1,2,3,4-tetrachlorodibenzodioxine REASON: CORRECT Matches '
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1cc2oc3c(Cl)c(Cl)c(Cl)c(Cl)c3c2c(Cl)c1Cl NAME: '
               '1,2,3,4,7,8,9-Heptachlorodibenzofuran REASON: CORRECT Matches '
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1cc2Oc3c(Cl)c(Cl)c(Cl)c(Cl)c3Oc2cc1Cl NAME: '
               '1,2,3,4,7,8-Hexachlorodibenzodioxin REASON: CORRECT Matches '
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1cc2oc3cc(Cl)c(Cl)c(Cl)c3c2c(Cl)c1Cl NAME: '
               '1,2,3,7,8,9-Hexachlorodibenzofuran REASON: CORRECT Matches '
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1cc2c(oc3c(Cl)c(Cl)c(Cl)c(Cl)c23)c(Cl)c1Cl NAME: '
               '1,2,3,4,6,7,8-Heptachlorodibenzofuran REASON: CORRECT Matches '
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1cc2Oc3c(Oc2c(Cl)c1Cl)cc(Cl)c(Cl)c3Cl NAME: '
               '1,2,3,6,7,8-Hexachlorodibenzodioxin REASON: CORRECT Matches '
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1cc(Cl)c(cc1Cl)-c1cc(Cl)c(Cl)cc1Cl NAME: '
               "2,2',4,4',5,5'-hexachlorobiphenyl REASON: CORRECT Matches "
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1cc2oc3c(Cl)c(Cl)c(Cl)c(Cl)c3c2cc1Cl NAME: '
               '1,2,3,4,7,8-Hexachlorodibenzofuran REASON: CORRECT Matches '
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1ccc(-c2cc(Cl)c(Cl)c(Cl)c2)c(Cl)c1Cl NAME: '
               "2,3,3',4,4',5'-Hexachlorobiphenyl REASON: CORRECT Matches "
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1ccc(cc1Cl)-c1cc(Cl)c(Cl)cc1Cl NAME: '
               "2,3',4,4',5-Pentachlorobiphenyl REASON: CORRECT Matches "
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1cc2oc3cc(Cl)c(Cl)c(Cl)c3c2cc1Cl NAME: '
               '1,2,3,7,8-Pentachlorodibenzofuran REASON: CORRECT Matches '
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: ClC=1C(C=2C(Cl)=CC(Cl)=C(Cl)C2)=CC(Cl)=C(Cl)C1Cl '
               'NAME: PCB180 REASON: CORRECT Matches polychlorinated '
               'dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1cc2Oc3cc(Cl)c(Cl)c(Cl)c3Oc2c(Cl)c1Cl NAME: '
               '1,2,3,7,8,9-Hexachlorodibenzodioxin REASON: CORRECT Matches '
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1cc2Oc3cc(Cl)c(Cl)cc3Oc2cc1Cl NAME: '
               '2,3,7,8-tetrachlorodibenzodioxine REASON: CORRECT Matches '
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1ccc(c(Cl)c1)-c1cc(Cl)c(Cl)c(Cl)c1 NAME: '
               "2,3',4,4',5'-Pentachlorobiphenyl REASON: CORRECT Matches "
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1ccc(cc1Cl)-c1ccc(Cl)c(Cl)c1Cl NAME: '
               "2,3,4,3',4'-Pentachlorobiphenyl REASON: CORRECT Matches "
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1c(Cl)c(Cl)c2Oc3c(Cl)c(Cl)c(Cl)c(Cl)c3Oc2c1Cl '
               'NAME: 1,2,3,4,6,7,8,9-Octachlorodibenzo-p-dioxin REASON: '
               'CORRECT Matches polychlorinated dibenzodioxin or related '
               'compound pattern.\n'
               ' * SMILES: ClC1=C(C=2C(Cl)=CC(Cl)=C(Cl)C2)C=CC(Cl)=C1Cl NAME: '
               'PCB138 REASON: CORRECT Matches polychlorinated dibenzodioxin '
               'or related compound pattern.\n'
               ' * SMILES: Clc1ccc(cc1Cl)-c1cc(Cl)c(Cl)c(Cl)c1Cl NAME: '
               "2,3,3',4,4',5-Hexachlorobiphenyl REASON: CORRECT Matches "
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1cc(Cl)c2Oc3cc(Cl)cc(Cl)c3Oc2c1 NAME: '
               '1,3,6,8-tetrachlorodibenzodioxine REASON: CORRECT Matches '
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1cc2Oc3cc(Cl)c(Cl)c(Cl)c3Oc2cc1Cl NAME: '
               '1,2,3,7,8-Pentachlorodibenzodioxin REASON: CORRECT Matches '
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               'False positives: SMILES: '
               'Oc1c(Cl)c(Cl)cc(c1Cl)-c1cc(Cl)c(Cl)c(Cl)c1Cl NAME: '
               "2,2',3',4,4',5,5'-Heptachloro-3-biphenylol REASON: WRONGLY "
               'CLASSIFIED Matches polychlorinated dibenzodioxin or related '
               'compound pattern.\n'
               ' * SMILES: Oc1c(Cl)cc2c(oc3cc(Cl)c(Cl)c(Cl)c23)c1Cl NAME: '
               '7-hydroxy-1,2,3,6,8-pentachlorodibenzofuran REASON: WRONGLY '
               'CLASSIFIED Matches polychlorinated dibenzodioxin or related '
               'compound pattern.\n'
               'False negatives: SMILES: Clc1ccc(cc1)-c1cc(Cl)c(Cl)c(Cl)c1Cl '
               "NAME: 2,3,4,4',5-Pentachlorobiphenyl REASON: MISSED Does not "
               'match polychlorinated dibenzodioxin or related compound '
               'pattern.\n'
               ' * SMILES: Clc1cccc(c1)-c1cc(Cl)c(Cl)c(Cl)c1Cl NAME: '
               "2,3,3',4,5-pentachlorobiphenyl REASON: MISSED Does not match "
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1cc(Cl)cc(c1)-c1cc(Cl)cc(Cl)c1 NAME: '
               "3,3',5,5'-tetrachlorobiphenyl REASON: MISSED Does not match "
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: '
               'ClC1=C(O)C=2C(=O)[C@]3(O)C(=O)C4=C(C5=C(OC)C=C(OC)C(=C5C)Cl)C=C(OC)C(=C4C[C@@H]3C(C2C=C1O)(C)C)Cl '
               'NAME: Formicamycin F REASON: MISSED Insufficient halogen atoms '
               'for a likely related compound.\n'
               ' * SMILES: Clc1c(Cl)c(Cl)c(-c2ccccc2)c(Cl)c1Cl NAME: '
               '2,3,4,5,6-pentachlorobiphenyl REASON: MISSED Does not match '
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1cc(Cl)cc(c1)-c1ccccc1 NAME: '
               '3,5-dichlorobiphenyl REASON: MISSED Insufficient halogen atoms '
               'for a likely related compound.\n'
               ' * SMILES: Clc1ccc(cc1)-c1c(Cl)cc(Cl)cc1Cl NAME: '
               "2,4,4',6-tetrachlorobiphenyl REASON: MISSED Does not match "
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1cccc(Cl)c1-c1ccccc1 NAME: 2,6-dichlorobiphenyl '
               'REASON: MISSED Insufficient halogen atoms for a likely related '
               'compound.\n'
               ' * SMILES: '
               'ClC1=C(O)C(Cl)=C2C([C@H]3CC=4C(Cl)=C(OC)C=C(C4C([C@]3(C(C2=C1O)=O)O)=O)C5=C(O)C=C(OC)C(=C5C)Cl)(C)C '
               'NAME: Formicamycin I REASON: MISSED Does not match '
               'polychlorinated dibenzodioxin or related compound pattern.\n'
               ' * SMILES: Clc1cccc(c1)-c1cc(Cl)ccc1Cl NAME: '
               "2,3',5-trichlorobiphenyl REASON: MISSED Insufficient halogen "
               'atoms for a likely related compound.\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1OC(=O)C=C1',
                                     'name': '(S)-5-methylfuran-2(5H)-one',
                                     'reason': 'Insufficient halogen atoms for '
                                               'a likely related compound.'},
                                 {   'smiles': 'Cc1cccc(c1)C(O)=O',
                                     'name': 'm-toluic acid',
                                     'reason': 'Insufficient halogen atoms for '
                                               'a likely related compound.'},
                                 {   'smiles': 'CN(C)CC(=O)NC[C@@H]1[C@@H]([C@@H](N1)CO)C2=CC=C(C=C2)C3=CC=CC(=C3)C#N',
                                     'name': 'N-[[(2S,3S,4R)-3-[4-(3-cyanophenyl)phenyl]-4-(hydroxymethyl)-2-azetidinyl]methyl]-2-(dimethylamino)acetamide',
                                     'reason': 'Insufficient halogen atoms for '
                                               'a likely related compound.'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OCCCCCCCCCC(=O)CC(O)=O)[C@H](O)C[C@H]1O',
                                     'name': 'bkos#20',
                                     'reason': 'Insufficient halogen atoms for '
                                               'a likely related compound.'},
                                 {   'smiles': 'C([C@@](COC(CCCCCCCCC/C=C\\CCCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])OP(=O)(O)OC[C@@](CO)([H])O',
                                     'name': 'PG(18:1(11Z)/18:1(11Z))',
                                     'reason': 'Insufficient halogen atoms for '
                                               'a likely related compound.'},
                                 {   'smiles': 'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC=2C3=C(C=CC=C3)NC2)CCC(=O)N)(C)C)(C)C)(CC)C)CCC1)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC4=CC=CC=C4)(C)C)CO)(C)C)(C)C)CC(C)C)CCC(=O)N)(C)C)(C)C)C)C)(C)C',
                                     'name': 'Chrysospermin B',
                                     'reason': 'Insufficient halogen atoms for '
                                               'a likely related compound.'},
                                 {   'smiles': 'CC1(C=CC2=C3C(=CC(=C2O1)OC)C(=O)C(=CO3)C4=CC5=C(C=C4OC)OCO5)C',
                                     'name': '6-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-8,8-dimethyl-4-pyrano[2,3-h][1]benzopyranone',
                                     'reason': 'Insufficient halogen atoms for '
                                               'a likely related compound.'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)NCC(O)=O)C',
                                     'name': 'Glu-Thr-Gly',
                                     'reason': 'Insufficient halogen atoms for '
                                               'a likely related compound.'},
                                 {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]3[C@@H](O)[C@H](O)[C@H](O[C@@H]3OC[C@H]4O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]5CO)O)[C@@H](O)[C@@H](O)[C@@H]4O)CO)[C@@H]6O[C@@H]([C@H](O)[C@H](O[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]6O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3S,4S,5S,6S)-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-2-(hydroxymethyl)-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'Insufficient halogen atoms for '
                                               'a likely related compound.'},
                                 {   'smiles': 'COCCOC(C)C(O)=O',
                                     'name': '2-(2-methoxyethoxy)propanoic '
                                             'acid',
                                     'reason': 'Insufficient halogen atoms for '
                                               'a likely related compound.'}],
    'sample_false_negatives': [   {   'smiles': 'Clc1cc(Cl)cc(c1)-c1cc(Cl)cc(Cl)c1',
                                      'name': "3,3',5,5'-tetrachlorobiphenyl",
                                      'reason': 'Does not match '
                                                'polychlorinated dibenzodioxin '
                                                'or related compound pattern.'},
                                  {   'smiles': 'Clc1cc2Oc3c(Cl)c(Cl)c(Cl)c(Cl)c3Oc2c(Cl)c1Cl',
                                      'name': '1,2,3,4,6,7,8-Heptachlorodibenzodioxin',
                                      'reason': 'Does not match '
                                                'polychlorinated dibenzodioxin '
                                                'or related compound pattern.'},
                                  {   'smiles': 'ClC1=C(O)C=2C(=O)[C@]3(O)C(=O)C4=C(C5=C(OC)C=C(OC)C(=C5C)Cl)C=C(OC)C(=C4C[C@@H]3C(C2C=C1O)(C)C)Cl',
                                      'name': 'Formicamycin F',
                                      'reason': 'Does not match '
                                                'polychlorinated dibenzodioxin '
                                                'or related compound pattern.'},
                                  {   'smiles': 'Clc1c(Cl)c(Cl)c2Oc3ccccc3Oc2c1Cl',
                                      'name': '1,2,3,4-tetrachlorodibenzodioxine',
                                      'reason': 'Does not match '
                                                'polychlorinated dibenzodioxin '
                                                'or related compound pattern.'},
                                  {   'smiles': 'Clc1cc(Cl)cc(c1)-c1ccccc1',
                                      'name': '3,5-dichlorobiphenyl',
                                      'reason': 'Insufficient halogen atoms '
                                                'for a likely related '
                                                'compound.'},
                                  {   'smiles': 'Clc1cccc(Cl)c1-c1ccccc1',
                                      'name': '2,6-dichlorobiphenyl',
                                      'reason': 'Insufficient halogen atoms '
                                                'for a likely related '
                                                'compound.'},
                                  {   'smiles': 'Clc1cc2Oc3c(Cl)c(Cl)c(Cl)c(Cl)c3Oc2cc1Cl',
                                      'name': '1,2,3,4,7,8-Hexachlorodibenzodioxin',
                                      'reason': 'Does not match '
                                                'polychlorinated dibenzodioxin '
                                                'or related compound pattern.'},
                                  {   'smiles': 'ClC1=C(O)C(Cl)=C2C([C@H]3CC=4C(Cl)=C(OC)C=C(C4C([C@]3(C(C2=C1O)=O)O)=O)C5=C(O)C=C(OC)C(=C5C)Cl)(C)C',
                                      'name': 'Formicamycin I',
                                      'reason': 'Does not match '
                                                'polychlorinated dibenzodioxin '
                                                'or related compound pattern.'},
                                  {   'smiles': 'OC1=C(C=C(C=C1[N+]([O-])=O)Cl)C2=C(O)C([N+](=O)[O-])=CC(=C2)Cl',
                                      'name': 'niclofolan',
                                      'reason': 'Insufficient halogen atoms '
                                                'for a likely related '
                                                'compound.'},
                                  {   'smiles': 'Clc1cc2Oc3c(Oc2c(Cl)c1Cl)cc(Cl)c(Cl)c3Cl',
                                      'name': '1,2,3,6,7,8-Hexachlorodibenzodioxin',
                                      'reason': 'Does not match '
                                                'polychlorinated dibenzodioxin '
                                                'or related compound '
                                                'pattern.'}],
    'attempt': 4,
    'reasoning': 'Based on the outcomes of the previous attempt, there are '
                 'several issues contributing to the low F1 score. Let’s '
                 'examine and address these:\n'
                 '\n'
                 '1. **False Positives**: The program misclassified molecules '
                 "as belonging to the target class when they shouldn't. The "
                 'presence of hydroxy groups (like -OH) seems to confuse the '
                 'pattern matching. We should ensure that SMARTS patterns '
                 'exclude molecules with functional groups that would '
                 'typically not belong in these types of compounds.\n'
                 '\n'
                 '2. **False Negatives**: The program failed to identify '
                 'several compounds that should belong to the target class. '
                 'This might indicate inadequate or overly specific SMARTS '
                 'patterns, or failure to encompass the variability within the '
                 'class (like insufficient bromine or chlorine atoms being too '
                 'strictly enforced).\n'
                 '\n'
                 '3. **SMARTS Pattern Flexibility**: The existing SMARTS '
                 'patterns are potentially too rigid. These compounds are '
                 'quite diverse, and patterns should be generalized to capture '
                 'this diversity without setting unrealistic constraints.\n'
                 '\n'
                 'To improve:\n'
                 '- Broaden the existing SMARTS patterns to catch the class '
                 'diversity.\n'
                 '- Add checks for high substitution of chlorine or bromine '
                 'atoms but relax the absolute threshold slightly to capture '
                 'more true positives.\n'
                 '- Combine specific substructural checks to exclude common '
                 'false positive confounders efficiently.\n'
                 '\n'
                 "Let's implement these improvements:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 28,
    'num_false_positives': 2,
    'num_true_negatives': 142242,
    'num_false_negatives': 28,
    'num_negatives': None,
    'precision': 0.9333333333333333,
    'recall': 0.5,
    'f1': 0.6511627906976745,
    'accuracy': 0.9997891777933943,
    'negative_predictive_value': 0.9998031911154847}