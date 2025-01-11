"""
Classifies: CHEBI:17334 penicillin
"""
"""
Classifies: CHEBI:17334 penicillin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    A penicillin is a member of the group of substituted penams containing:
    - Two methyl substituents at position 2
    - A carboxylate substituent at position 3
    - A carboxamido group at position 6

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penicillin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the penicillin core structure pattern
    penicillin_core = Chem.MolFromSmarts("[C@@H]1[C@H]2N(C1=O)[C@H](C(S2)(C)C)C(=O)O")
    if not mol.HasSubstructMatch(penicillin_core):
        return False, "No penicillin core structure found"

    # Check for two methyl groups at position 2
    methyl_pattern = Chem.MolFromSmarts("[C@@H]1[C@H]2N(C1=O)[C@H](C(S2)(C)C)C(=O)O")
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)
    if len(methyl_matches) == 0:
        return False, "No methyl groups found at position 2"

    # Check for carboxylate group at position 3
    carboxylate_pattern = Chem.MolFromSmarts("[C@@H]1[C@H]2N(C1=O)[C@H](C(S2)(C)C)C(=O)O")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(carboxylate_matches) == 0:
        return False, "No carboxylate group found at position 3"

    # Check for carboxamido group at position 6
    carboxamido_pattern = Chem.MolFromSmarts("[C@@H]1[C@H]2N(C1=O)[C@H](C(S2)(C)C)C(=O)O")
    carboxamido_matches = mol.GetSubstructMatches(carboxamido_pattern)
    if len(carboxamido_matches) == 0:
        return False, "No carboxamido group found at position 6"

    return True, "Contains penicillin core structure with required substituents"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17334',
                          'name': 'penicillin',
                          'definition': 'Any member of the group of substituted penams containing two methyl substituents at position 2, a carboxylate substituent at position 3 and a carboxamido group at position 6.',
                          'parents': ['CHEBI:17334', 'CHEBI:17334']},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17334',
                          'name': 'penicillin',
                          'definition': 'Any member of the group of '
                                        'substituted penams containing two '
                                        'methyl substituents at position 2, a '
                                        'carboxylate substituent at position 3 '
                                        'and a carboxamido group at position '
                                        '6.',
                          'parents': ['CHEBI:25865'],
                          'xrefs': [   'KEGG:C00395',
                                       'PMID:11851248',
                                       'PMID:12833570',
                                       'PMID:1502708',
                                       'PMID:16033609',
                                       'PMID:7061385',
                                       'PMID:7798534',
                                       'Wikipedia:Penicillin'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': '[H]\\C(C)=C1/CN(C)[C@]2([H])Cc3c(CC[C@]1([H])[C@]2([H])C)[nH]c1ccccc31',
                                     'name': 'vobasan',
                                     'reason': 'No penicillin core structure '
                                               'found'},
                                 {   'smiles': 'CC(O)C(O)=O',
                                     'name': '2-hydroxypropanoic acid',
                                     'reason': 'No penicillin core structure '
                                               'found'},
                                 {   'smiles': 'O=C(O)C[C@@](O)(CC/C=C(/COC(=O)C)\\C)C',
                                     'name': 'Penicimonoterpene',
                                     'reason': 'No penicillin core structure '
                                               'found'},
                                 {   'smiles': 'O=C1NC(O)(CCCC)[C@@H]2[C@H]1O2',
                                     'name': '(+)-epogymnolactam',
                                     'reason': 'No penicillin core structure '
                                               'found'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCC(=O)C[C@@H]1[C@H]([C@H](O)C[C@@H]1O)/C=C/[C@@H](O)CCCCC)COC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)(O)(O)=O',
                                     'name': 'PA(22:6(4Z,7Z,10Z,13Z,16Z,19Z)/6 '
                                             'keto-PGF1alpha)',
                                     'reason': 'No penicillin core structure '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](O[C@H]2CC[C@@]3(C)[C@@H](CC[C@]4(C)[C@@H]3CC=C3[C@@H]5CC(C)(C)[C@@H](O)[C@H](O)[C@]5(CO)C[C@H](O)[C@@]43C)C2(C)C)[C@H](O)[C@H](O)[C@H]1O[C@@H]1OC[C@H](O)[C@H](O)[C@H]1O',
                                     'name': 'Pittoside A',
                                     'reason': 'No penicillin core structure '
                                               'found'},
                                 {   'smiles': 'C1CC(C1)C(=O)NC2=CC3=C(C=C2)O[C@@H]4[C@H]3C[C@@H](O[C@H]4CO)CC(=O)NCCC5=CC=CC=C5',
                                     'name': 'N-[(1S,3R,4aS,9aR)-1-(hydroxymethyl)-3-[2-oxo-2-(2-phenylethylamino)ethyl]-3,4,4a,9a-tetrahydro-1H-pyrano[3,4-b]benzofuran-6-yl]cyclobutanecarboxamide',
                                     'reason': 'No penicillin core structure '
                                               'found'},
                                 {   'smiles': 'O1C(C(C2=C1C(OC)=CC(=C2)/C=C/CO)C)C3=CC(OC)=C(OC)C=C3',
                                     'name': 'Fragransol C',
                                     'reason': 'No penicillin core structure '
                                               'found'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O',
                                     'name': 'PS(20:0/17:0)',
                                     'reason': 'No penicillin core structure '
                                               'found'},
                                 {   'smiles': 'C[C@H]1CN(C(=O)C2=C(C(=CC=C2)NC(=O)CCCCCC(=O)NC3=CC=CC=C3N)O[C@H]1CN(C)S(=O)(=O)C4=CC=C(C=C4)Cl)[C@@H](C)CO',
                                     'name': "N-(2-aminophenyl)-N'-[(2R,3S)-2-[[(4-chlorophenyl)sulfonyl-methylamino]methyl]-5-[(2S)-1-hydroxypropan-2-yl]-3-methyl-6-oxo-3,4-dihydro-2H-1,5-benzoxazocin-10-yl]heptanediamide",
                                     'reason': 'No penicillin core structure '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@]2(NC(=O)C([H])(C(O)=O)c1ccsc1)OC)C(O)=O',
                                      'name': 'temocillin',
                                      'reason': 'No penicillin core structure '
                                                'found'}],
    'attempt': 0,
    'reasoning': 'To classify a molecule as a penicillin, we need to identify '
                 'the core structure of penicillin, which includes a bicyclic '
                 'system (4-thia-1-azabicyclo[3.2.0]heptane) with specific '
                 'substituents: two methyl groups at position 2, a carboxylate '
                 'group at position 3, and a carboxamido group at position 6. '
                 'The function `is_penicillin` will check for these structural '
                 'features using RDKit.\n'
                 '\n'
                 'Here’s the Python code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 49,
    'num_false_positives': 75,
    'num_true_negatives': 142175,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.3951612903225806,
    'recall': 0.98,
    'f1': 0.5632183908045977,
    'accuracy': 0.9994659170765987,
    'negative_predictive_value': 0.9999929664641009}