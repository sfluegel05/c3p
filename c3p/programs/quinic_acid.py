"""
Classifies: CHEBI:26493 quinic acid
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cyclohexane core with multiple hydroxyls and one carboxylic acid
    # Quinic acid base structure: OC1CC(O)(C(CCC1O)C(O)=O)O
    quinic_acid_base_pattern = Chem.MolFromSmarts("OC1CCC(O)(C(O)=O)CC1O")
    
    if not mol.HasSubstructMatch(quinic_acid_base_pattern):
        return False, "Mismatch in quinic acid core structure (cyclohexane with hydroxyls and carboxylic acid)"

    # If matches, check for any additional features or derivatives
    # Note: This part is flexible as quinic acid derivatives can vary
    # For simple quinic acids
    base_match = mol.HasSubstructMatch(quinic_acid_base_pattern)
    if base_match:
        return True, "Matches quinic acid core structure"

    return False, "Does not match essential quinic acid features"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26493',
                          'name': 'quinic acid',
                          'definition': 'A cyclitol carboxylic acid.',
                          'parents': ['CHEBI:36123'],
                          'xrefs': [   'PMID:15589483',
                                       'PMID:19952409',
                                       'PMID:23590498'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'Mismatch in quinic acid core '
                                               'structure (cyclohexane with '
                                               'hydroxyls and carboxylic '
                                               'acid)'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'Mismatch in quinic acid core '
                                               'structure (cyclohexane with '
                                               'hydroxyls and carboxylic '
                                               'acid)'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'Mismatch in quinic acid core '
                                               'structure (cyclohexane with '
                                               'hydroxyls and carboxylic '
                                               'acid)'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'Mismatch in quinic acid core '
                                               'structure (cyclohexane with '
                                               'hydroxyls and carboxylic '
                                               'acid)'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'Mismatch in quinic acid core '
                                               'structure (cyclohexane with '
                                               'hydroxyls and carboxylic '
                                               'acid)'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': 'Mismatch in quinic acid core '
                                               'structure (cyclohexane with '
                                               'hydroxyls and carboxylic '
                                               'acid)'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'Mismatch in quinic acid core '
                                               'structure (cyclohexane with '
                                               'hydroxyls and carboxylic '
                                               'acid)'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'Mismatch in quinic acid core '
                                               'structure (cyclohexane with '
                                               'hydroxyls and carboxylic '
                                               'acid)'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': 'Mismatch in quinic acid core '
                                               'structure (cyclohexane with '
                                               'hydroxyls and carboxylic '
                                               'acid)'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'Mismatch in quinic acid core '
                                               'structure (cyclohexane with '
                                               'hydroxyls and carboxylic '
                                               'acid)'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To determine if a molecule is a quinic acid based on its '
                 "SMILES string, we'll focus on identifying key structural "
                 'features from the provided examples. Quinic acids are '
                 'characterized by a cyclohexane ring with multiple hydroxyl '
                 'groups and a carboxylic acid group. Often, this structure is '
                 'involved in ester or conjugated bonds in various '
                 'derivatives.\n'
                 '\n'
                 "Here are the key features we'll look for to classify a "
                 'compound as a quinic acid:\n'
                 '\n'
                 '1. A cyclohexane ring with hydroxyl groups.\n'
                 '2. A carboxylic acid group attached to the cyclohexane.\n'
                 '3. The possibility of ester linkages or other functional '
                 'groups attached to the hydroxyl positions.\n'
                 '\n'
                 'We can use SMARTS patterns to search for these features in '
                 'the molecule.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 25,
    'num_false_positives': 32,
    'num_true_negatives': 142243,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.43859649122807015,
    'recall': 1.0,
    'f1': 0.6097560975609756,
    'accuracy': 0.9997751229796206,
    'negative_predictive_value': 1.0}