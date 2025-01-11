"""
Classifies: CHEBI:28863 flavanones
"""
from rdkit import Chem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone based on its SMILES string.
    A flavanone has a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a generalized SMARTS pattern for the core flavanone skeleton
    # Focus on a 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one structure
    flavanone_pattern = Chem.MolFromSmarts("c1ccccc1C2CC(=O)c3ccccc3O2")
    # Core feature: aromatic benzene ring, carbonyl on chromanone, 2nd aryl group attached on chromanone
    
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone core skeleton found"
    
    # Additional checks for variability in substituents can be done if needed
    
    return True, "Contains 3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:28863',
                          'name': 'flavanones',
                          'definition': 'Members of the class of flavans with '
                                        'a '
                                        '3,4-dihydro-2-aryl-2H-1-benzopyran-4-one '
                                        'skeleton and its substituted '
                                        'derivatives.',
                          'parents': ['CHEBI:38672', 'CHEBI:3992'],
                          'xrefs': [   'MetaCyc:FLAVANONES',
                                       'Wikipedia:Flavanone'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               "Error: name 'flayanone_pattern' is not defined\n"
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
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1OC(=O)C=C1',
                                     'name': '(S)-5-methylfuran-2(5H)-one',
                                     'reason': 'No flavanone core skeleton '
                                               'found'},
                                 {   'smiles': 'Cc1cccc(c1)C(O)=O',
                                     'name': 'm-toluic acid',
                                     'reason': 'No flavanone core skeleton '
                                               'found'},
                                 {   'smiles': 'CN(C)CC(=O)NC[C@@H]1[C@@H]([C@@H](N1)CO)C2=CC=C(C=C2)C3=CC=CC(=C3)C#N',
                                     'name': 'N-[[(2S,3S,4R)-3-[4-(3-cyanophenyl)phenyl]-4-(hydroxymethyl)-2-azetidinyl]methyl]-2-(dimethylamino)acetamide',
                                     'reason': 'No flavanone core skeleton '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OCCCCCCCCCC(=O)CC(O)=O)[C@H](O)C[C@H]1O',
                                     'name': 'bkos#20',
                                     'reason': 'No flavanone core skeleton '
                                               'found'},
                                 {   'smiles': 'C([C@@](COC(CCCCCCCCC/C=C\\CCCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])OP(=O)(O)OC[C@@](CO)([H])O',
                                     'name': 'PG(18:1(11Z)/18:1(11Z))',
                                     'reason': 'No flavanone core skeleton '
                                               'found'},
                                 {   'smiles': 'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC=2C3=C(C=CC=C3)NC2)CCC(=O)N)(C)C)(C)C)(CC)C)CCC1)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC4=CC=CC=C4)(C)C)CO)(C)C)(C)C)CC(C)C)CCC(=O)N)(C)C)(C)C)C)C)(C)C',
                                     'name': 'Chrysospermin B',
                                     'reason': 'No flavanone core skeleton '
                                               'found'},
                                 {   'smiles': 'CC1(C=CC2=C3C(=CC(=C2O1)OC)C(=O)C(=CO3)C4=CC5=C(C=C4OC)OCO5)C',
                                     'name': '6-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-8,8-dimethyl-4-pyrano[2,3-h][1]benzopyranone',
                                     'reason': 'No flavanone core skeleton '
                                               'found'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)NCC(O)=O)C',
                                     'name': 'Glu-Thr-Gly',
                                     'reason': 'No flavanone core skeleton '
                                               'found'},
                                 {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]3[C@@H](O)[C@H](O)[C@H](O[C@@H]3OC[C@H]4O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]5CO)O)[C@@H](O)[C@@H](O)[C@@H]4O)CO)[C@@H]6O[C@@H]([C@H](O)[C@H](O[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]6O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3S,4S,5S,6S)-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-2-(hydroxymethyl)-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No flavanone core skeleton '
                                               'found'},
                                 {   'smiles': 'COCCOC(C)C(O)=O',
                                     'name': '2-(2-methoxyethoxy)propanoic '
                                             'acid',
                                     'reason': 'No flavanone core skeleton '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'COc1ccc(C2COc3cc(O)cc(O)c3C2=O)c(OC)c1',
                                      'name': 'homoferreirin',
                                      'reason': 'No flavanone core skeleton '
                                                'found'},
                                  {   'smiles': 'COC1=CC=C(C=C1)C1=CC(=O)C2=C(O)C(OC)=C(O[C@@H]3O[C@H](CO[C@@H]4O[C@@H](C)[C@H](O)[C@@H](O)[C@H]4O)[C@@H](O)[C@H](O)[C@H]3O)C=C2O1',
                                      'name': 'pectolinarin',
                                      'reason': 'No flavanone core skeleton '
                                                'found'},
                                  {   'smiles': 'COc1ccc(cc1O)-c1cc(=O)c2c(O)cc(O[C@@H]3O[C@H](CO[C@@H]4O[C@@H](C)[C@H](O)[C@@H](O)[C@H]4O)[C@@H](O)[C@H](O)[C@H]3O)cc2o1',
                                      'name': 'diosmin',
                                      'reason': 'No flavanone core skeleton '
                                                'found'}],
    'attempt': 5,
    'reasoning': 'Based on your previous code and its outcomes, it appears '
                 'that there was a simple typo in the variable name '
                 '`flayanone_pattern`, which should be `flavanone_pattern`. '
                 'Since this essential line contained the error, the function '
                 'was unable to run properly, as the correct substructure '
                 "search could not be executed. Furthermore, let's review if "
                 'the SMARTS pattern accurately captures the descriptor of a '
                 'flavanone, ensuring the pattern covers the core flavanone '
                 'structure.\n'
                 '\n'
                 'The fundamental structure of flavanones involves a '
                 '3,4-dihydro-2-aryl-2H-1-benzopyran-4-one framework. This can '
                 'be represented by focusing on:\n'
                 '- A benzene ring\n'
                 '- A central carbon in a dihydrochromanone (this accounts for '
                 'the 3,4-dihydro-2H-pyran component)\n'
                 '- Carbonyl group at the 4-position of the chromanone\n'
                 '- Another aryl (phenyl) group attached\n'
                 '\n'
                 "Let's create a function with these considerations, using the "
                 'corrected pattern name:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 279,
    'num_false_positives': 100,
    'num_true_negatives': 55253,
    'num_false_negatives': 3,
    'num_negatives': None,
    'precision': 0.7361477572559367,
    'recall': 0.9893617021276596,
    'f1': 0.844175491679274,
    'accuracy': 0.9981486474341691,
    'negative_predictive_value': 0.999945707253511}