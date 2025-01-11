"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for ketone group (C=O) at the second carbon
    ketone_pattern = Chem.MolFromSmarts("CC(=O)")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "Missing 2-oxo group (ketone group at the second carbon)"
        
    # Look for carboxylic acid group (C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Missing carboxylic acid group"

    return True, "Contains 2-oxo group and carboxylic acid group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35910',
                          'name': '2-oxo monocarboxylic acid',
                          'definition': 'Any monocarboxylic acid having a '
                                        '2-oxo substituent.',
                          'parents': ['CHEBI:35871'],
                          'xrefs': ['KEGG:C00161'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'Missing carboxylic acid group'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'Missing carboxylic acid group'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'Missing 2-oxo group (ketone '
                                               'group at the second carbon)'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'Missing carboxylic acid group'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'Missing 2-oxo group (ketone '
                                               'group at the second carbon)'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)NC3=CC=C(C=C3)Cl)C(=O)N(C[C@@H]1OC)C)C)C',
                                     'name': '1-(4-chlorophenyl)-3-[(4S,7R,8R)-8-methoxy-4,5,7,10-tetramethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]urea',
                                     'reason': 'Missing 2-oxo group (ketone '
                                               'group at the second carbon)'},
                                 {   'smiles': 'C1C[C@@H]([C@@H](O[C@@H]1CCNC(=O)C2CCOCC2)CO)NC(=O)NC3=CC(=CC=C3)Cl',
                                     'name': 'N-[2-[(2S,5S,6R)-5-[[(3-chloroanilino)-oxomethyl]amino]-6-(hydroxymethyl)-2-oxanyl]ethyl]-4-oxanecarboxamide',
                                     'reason': 'Missing carboxylic acid group'},
                                 {   'smiles': '[H][C@@]1(CNc2nc(N)[nH]c(=O)c2N1)[C@@H](O)[C@H](C)O',
                                     'name': 'sapropterin',
                                     'reason': 'Missing 2-oxo group (ketone '
                                               'group at the second carbon)'},
                                 {   'smiles': 'CCCCN(CCCC)CCNC(=O)CNC(=O)C1=CC=CS1',
                                     'name': 'N-[2-[2-(dibutylamino)ethylamino]-2-oxoethyl]-2-thiophenecarboxamide',
                                     'reason': 'Missing carboxylic acid group'},
                                 {   'smiles': 'C[C@@H](C1=CC=CC=C1)NC(=O)C[C@@H]2CC[C@H]([C@H](O2)CO)NC(=O)C3CCC3',
                                     'name': 'N-[(2S,3R,6S)-2-(hydroxymethyl)-6-[2-oxo-2-[[(1S)-1-phenylethyl]amino]ethyl]-3-oxanyl]cyclobutanecarboxamide',
                                     'reason': 'Missing carboxylic acid '
                                               'group'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To classify a molecule as a 2-oxo monocarboxylic acid, we '
                 'need to identify two key structural features from its SMILES '
                 'string: \n'
                 '\n'
                 '1. The presence of a 2-oxo group, which typically represents '
                 'a ketone group positioned at the second carbon atom of the '
                 'molecule. This would appear in SMILES as "C(=O)", but '
                 'specifically must be at the second carbon in the chain.\n'
                 '2. The presence of a carboxylic acid group, represented as '
                 '"C(=O)O".\n'
                 '\n'
                 "Given that the description provided is 'None' for this "
                 'class, our task will be to match these functional groups to '
                 'identify whether the input SMILES string represents a 2-oxo '
                 'monocarboxylic acid. We will use the RDKit library to parse '
                 'and analyze the SMILES string. The function will return True '
                 'if both criteria are met, otherwise False.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 92,
    'num_false_positives': 100,
    'num_true_negatives': 164,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.4791666666666667,
    'recall': 1.0,
    'f1': 0.647887323943662,
    'accuracy': 0.7191011235955056,
    'negative_predictive_value': 1.0}