"""
Classifies: CHEBI:24043 flavones
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_flavones(smiles: str) -> (bool, str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    A flavone has a 2-aryl-1-benzopyran-4-one skeleton and its substituted derivatives.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more flexible core structure for flavones: 2-aryl-1-benzopyran-4-one
    flavone_core_pattern = Chem.MolFromSmarts("c1cc(ccc1)-c2cc(=O)c3ccccc3o2")
    
    # Ensure pattern was defined correctly
    if flavone_core_pattern is None:
        return None, "Failed to construct flavone core pattern"

    # Ensure molecule has flavone core structure
    if not mol.HasSubstructMatch(flavone_core_pattern):
        return False, "No flavone core structure found"

    return True, "Contains flavone core structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24043',
                          'name': 'flavones',
                          'definition': 'A member of the class of flavonoid '
                                        'with a 2-aryl-1-benzopyran-4-one '
                                        '(2-arylchromen-4-one) skeleton and '
                                        'its substituted derivatives.',
                          'parents': ['CHEBI:192499', 'CHEBI:47916'],
                          'xrefs': ['MetaCyc:Flavones', 'Wikipedia:Flavone'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'OC[C@H]1O[C@H]([C@H](O)[C@@H](O)[C@@H]1O)c1c(O)c([C@@H]2OC[C@H](O)[C@H](O)[C@H]2O)c(O)c2c1oc(cc2=O)-c1ccc(O)cc1 '
               'NAME: isoschaftoside REASON: MISSED No flavone core structure '
               'found\n'
               ' * SMILES: '
               'COc1cc(cc(OC)c1OC)-c1cc(=O)c2c(OC)c(OC)c(OC)c(OC)c2o1 NAME: '
               '2-(3,4,5-trimethoxyphenyl)-5,6,7,8-tetramethoxy-4H-1-benzopyran-4-one '
               'REASON: MISSED No flavone core structure found\n'
               ' * SMILES: '
               'COc1cc(ccc1O)-c1oc2cc(O[C@H]3O[C@@H](CO)[C@H](O)[C@@H](O)[C@@H]3O)cc(O)c2c(=O)c1O '
               'NAME: isorhamnetin 7-O-beta-L-glucopyranoside REASON: MISSED '
               'No flavone core structure found\n'
               ' * SMILES: O1C=2C(=C(O)C(=C(O)C2C=O)C)C(=O)C=C1C3=CC=CC=C3 '
               'NAME: Isounonal REASON: MISSED No flavone core structure '
               'found\n'
               ' * SMILES: '
               'OC[C@H]1OC(Oc2c(oc3ccccc3c2=O)-c2ccccc2)[C@H](O)[C@@H](O)[C@H]1O '
               'NAME: flavonol 3-O-D-galactoside REASON: MISSED No flavone '
               'core structure found\n'
               ' * SMILES: '
               'O[C@@H]1[C@@H](COC(=O)\\C=C\\c2ccc(O)cc2)O[C@@H](Oc2c(oc3cc(O)cc(O)c3c2=O)-c2ccc(O)cc2)[C@H](O)[C@H]1O '
               'NAME: tribuloside REASON: MISSED No flavone core structure '
               'found\n'
               ' * SMILES: '
               'O1C2=C(C(O)=C(CC=C(C)C)C(O)=C2)C(=O)C(O)=C1C3=CC=C(O)C=C3 '
               'NAME: Licoflavonol REASON: MISSED No flavone core structure '
               'found\n'
               ' * SMILES: '
               'OC[C@H]1O[C@@H](O[C@H]2[C@H](O[C@H]3[C@H](Oc4c(oc5cc(O)cc(O)c5c4=O)-c4ccc(O)cc4)O[C@H](CO)[C@@H](O)[C@@H]3O)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@@H]1O '
               'NAME: kaempferol '
               '3-O-beta-D-glucosyl-(1->2)-beta-D-glucosyl-(1->2)-beta-D-glucoside '
               'REASON: MISSED No flavone core structure found\n'
               ' * SMILES: ClC1=C(O)C2=C(OC(C3=CC=C(O)C=C3)=CC2=O)C=C1O NAME: '
               '6-Chloroapigenin REASON: MISSED No flavone core structure '
               'found\n'
               ' * SMILES: '
               'OC[C@H]1O[C@H]([C@H](O)[C@@H](O)[C@@H]1O)c1c(O)c([C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)c2oc(cc(=O)c2c1O)-c1ccc(O)c(O)c1 '
               'NAME: lucenin-2 REASON: MISSED No flavone core structure '
               'found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1OC(=O)C=C1',
                                     'name': '(S)-5-methylfuran-2(5H)-one',
                                     'reason': 'No flavone core structure '
                                               'found'},
                                 {   'smiles': 'Cc1cccc(c1)C(O)=O',
                                     'name': 'm-toluic acid',
                                     'reason': 'No flavone core structure '
                                               'found'},
                                 {   'smiles': 'CN(C)CC(=O)NC[C@@H]1[C@@H]([C@@H](N1)CO)C2=CC=C(C=C2)C3=CC=CC(=C3)C#N',
                                     'name': 'N-[[(2S,3S,4R)-3-[4-(3-cyanophenyl)phenyl]-4-(hydroxymethyl)-2-azetidinyl]methyl]-2-(dimethylamino)acetamide',
                                     'reason': 'No flavone core structure '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OCCCCCCCCCC(=O)CC(O)=O)[C@H](O)C[C@H]1O',
                                     'name': 'bkos#20',
                                     'reason': 'No flavone core structure '
                                               'found'},
                                 {   'smiles': 'C([C@@](COC(CCCCCCCCC/C=C\\CCCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])OP(=O)(O)OC[C@@](CO)([H])O',
                                     'name': 'PG(18:1(11Z)/18:1(11Z))',
                                     'reason': 'No flavone core structure '
                                               'found'},
                                 {   'smiles': 'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC=2C3=C(C=CC=C3)NC2)CCC(=O)N)(C)C)(C)C)(CC)C)CCC1)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC4=CC=CC=C4)(C)C)CO)(C)C)(C)C)CC(C)C)CCC(=O)N)(C)C)(C)C)C)C)(C)C',
                                     'name': 'Chrysospermin B',
                                     'reason': 'No flavone core structure '
                                               'found'},
                                 {   'smiles': 'CC1(C=CC2=C3C(=CC(=C2O1)OC)C(=O)C(=CO3)C4=CC5=C(C=C4OC)OCO5)C',
                                     'name': '6-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-8,8-dimethyl-4-pyrano[2,3-h][1]benzopyranone',
                                     'reason': 'No flavone core structure '
                                               'found'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)NCC(O)=O)C',
                                     'name': 'Glu-Thr-Gly',
                                     'reason': 'No flavone core structure '
                                               'found'},
                                 {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]3[C@@H](O)[C@H](O)[C@H](O[C@@H]3OC[C@H]4O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]5CO)O)[C@@H](O)[C@@H](O)[C@@H]4O)CO)[C@@H]6O[C@@H]([C@H](O)[C@H](O[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]6O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3S,4S,5S,6S)-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-2-(hydroxymethyl)-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No flavone core structure '
                                               'found'},
                                 {   'smiles': 'COCCOC(C)C(O)=O',
                                     'name': '2-(2-methoxyethoxy)propanoic '
                                             'acid',
                                     'reason': 'No flavone core structure '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'C1(OC2=CC(O)=C3C(OC=C(C3=O)C4=CC=C(O)C=C4)=C2)OC(COC(C)=O)C(O)C(O)C1O',
                                      'name': "6''-O-acetylgenistin",
                                      'reason': 'No flavone core structure '
                                                'found'},
                                  {   'smiles': '[C@H]1(C=2C3=C(C(C(C4=CC=C(O)C=C4)=CO3)=O)C=CC2[O-])[C@@H]([C@H]([C@H](O)[C@@H](CO)O1)O)O',
                                      'name': 'puerarin(1-)',
                                      'reason': 'No flavone core structure '
                                                'found'},
                                  {   'smiles': 'OC[C@@]1(O)CO[C@@H](OC[C@H]2O[C@H]([C@H](O)[C@@H](O)[C@@H]2O)c2c(O)ccc3c2occ(-c2ccc(O)cc2)c3=O)[C@@H]1O',
                                      'name': 'mirificin',
                                      'reason': 'No flavone core structure '
                                                'found'},
                                  {   'smiles': 'OC[C@H]1O[C@@H](Oc2cc(O)c3C(=O)C[C@H](Oc3c2[C@@H]2OC[C@H](O)[C@H](O)[C@H]2O)c2ccc(O)cc2)[C@H](O)[C@@H](O)[C@@H]1O',
                                      'name': '8-C-alpha-L-arabinopyranosyl-7-O-beta-D-glucopyranosylnaringenin',
                                      'reason': 'No flavone core structure '
                                                'found'},
                                  {   'smiles': 'Oc1ccc(cc1)[C@H]1Oc2cc(O)cc(O)c2C(=O)[C@@H]1[C@H]1[C@@H](Oc2cc(O)cc(O)c2C1=O)c1ccc(O)cc1',
                                      'name': 'isochamaejasmin',
                                      'reason': 'No flavone core structure '
                                                'found'},
                                  {   'smiles': 'COC1=CC2OC(=CC(=O)C2C(O)=C1[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)c1ccc(O)c(O)c1',
                                      'name': 'swertiajaponin',
                                      'reason': 'No flavone core structure '
                                                'found'},
                                  {   'smiles': 'C[C@@H]1O[C@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2Oc2ccc3C(=O)CC(Oc3c2)c2ccccc2)[C@H](O)[C@H](O)[C@H]1O',
                                      'name': '7-hydroxyflavanone '
                                              '7-O-[beta-L-rhamnosyl-(1->2)-beta-D-glucoside]',
                                      'reason': 'No flavone core structure '
                                                'found'},
                                  {   'smiles': 'O[C@@H]1[C@@H](O)[C@@H](O[C@@H]([C@H]1O)C(O)=O)Oc1cc(O)c2C(=O)C[C@@H](Oc2c1)c1ccc(O)c(O)c1',
                                      'name': '(2R)-eriodictoyl-7-O-beta-D-glucopyranosiduronic '
                                              'acid',
                                      'reason': 'No flavone core structure '
                                                'found'}],
    'attempt': 4,
    'reasoning': 'The previous attempt to classify flavones using SMILES '
                 'strings failed mainly due to the incorrect or too simplistic '
                 'SMARTS pattern used to identify the flavone core structure. '
                 'Flavones are characterized by having a '
                 '2-aryl-1-benzopyran-4-one skeleton, but this may include '
                 'various substitutions, and the previous pattern may not have '
                 'been comprehensive enough to capture all these variations '
                 'within the flavone class.\n'
                 '\n'
                 'To improve the classification of flavones, we should:\n'
                 '1. Use a more flexible and descriptive SMARTS pattern that '
                 'captures the essence of the flavone skeleton while allowing '
                 'for common variations and substitutions.\n'
                 '2. Consider the possibility of additional ring structures or '
                 'substitutions that may occur on flavones.\n'
                 '\n'
                 "Now, let's attempt to improve the function by defining a "
                 "better SMARTS pattern for flavones. Here's the modified "
                 'version:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 509,
    'num_false_positives': 100,
    'num_true_negatives': 13292,
    'num_false_negatives': 8,
    'num_negatives': None,
    'precision': 0.8357963875205254,
    'recall': 0.9845261121856866,
    'f1': 0.9040852575488456,
    'accuracy': 0.9922352433676037,
    'negative_predictive_value': 0.9993984962406015}