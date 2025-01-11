"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid contains two carboxylic acid groups (-COOH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dicarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group pattern (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    
    # Count carboxylic acid groups
    if len(carboxylic_acid_matches) < 2:
        return False, f"Found {len(carboxylic_acid_matches)} carboxylic acid groups, need at least 2"

    return True, "Contains at least two carboxylic acid groups"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35692',
                          'name': 'dicarboxylic acid',
                          'definition': 'Any carboxylic acid containing two '
                                        'carboxy groups.',
                          'parents': ['CHEBI:131927', 'CHEBI:33575'],
                          'xrefs': ['KEGG:C02028'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'Found 1 carboxylic acid '
                                               'groups, need at least 2'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'Found 0 carboxylic acid '
                                               'groups, need at least 2'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'Found 0 carboxylic acid '
                                               'groups, need at least 2'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'Found 0 carboxylic acid '
                                               'groups, need at least 2'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'Found 0 carboxylic acid '
                                               'groups, need at least 2'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'Found 1 carboxylic acid '
                                               'groups, need at least 2'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': 'Found 1 carboxylic acid '
                                               'groups, need at least 2'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'Found 0 carboxylic acid '
                                               'groups, need at least 2'},
                                 {   'smiles': 'SC[C@H](N)C(=O)N[C@H](C(=O)N[C@@H](CCC(=O)N)C(O)=O)CO',
                                     'name': 'Cys-Ser-Gln',
                                     'reason': 'Found 1 carboxylic acid '
                                               'groups, need at least 2'},
                                 {   'smiles': 'S(C(C)C(O)=O)CC1=CC=CC=C1',
                                     'name': '2-(Benzylthio)propanoic acid',
                                     'reason': 'Found 1 carboxylic acid '
                                               'groups, need at least 2'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To classify a molecule as a dicarboxylic acid using its '
                 'SMILES string, we need to identify the presence of two '
                 'carboxyl groups (-COOH). A dicarboxylic acid is '
                 'characterized by having at least two carboxylic acid groups '
                 'within its structure.\n'
                 '\n'
                 "Let's devise a function `is_dicarboxylic_acid` that will "
                 'take a SMILES string as input and determine if the molecule '
                 'represented by the string matches this chemical class. The '
                 'function should:\n'
                 '\n'
                 '1. Parse the SMILES string to create an RDKit molecule '
                 'object.\n'
                 '2. Identify and count the carboxylic acid groups within the '
                 'molecule using a SMARTS pattern.\n'
                 '3. Ensure that there are at least two such groups.\n'
                 '4. Return the classification as a boolean value along with a '
                 'reason for the classification.\n'
                 '\n'
                 'Below is the implementation of the `is_dicarboxylic_acid` '
                 'function.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 583,
    'num_false_positives': 100,
    'num_true_negatives': 536,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.8535871156661786,
    'recall': 1.0,
    'f1': 0.9210110584518167,
    'accuracy': 0.9179655455291222,
    'negative_predictive_value': 1.0}