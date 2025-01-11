"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid is defined as a fatty acid carrying one or more hydroxy substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O,H]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for hydroxy group (-OH)
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No hydroxy groups found"
    
    return True, "Contains carboxylic acid group with one or more hydroxy substituents"

__metadata__ = {
    'chemical_class': {
        'name': 'hydroxy fatty acid',
        'definition': 'Any fatty acid carrying one or more hydroxy substituents.'
    }
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24654',
                          'name': 'hydroxy fatty acid',
                          'definition': 'Any fatty acid carrying one or more '
                                        'hydroxy substituents.',
                          'parents': ['CHEBI:35366', 'CHEBI:35868'],
                          'xrefs': [   'LIPID_MAPS_class:LMFA0105',
                                       'PMID:18296335',
                                       'PMID:6419288',
                                       'PMID:8274032'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O=C1C2=C(OC(=C1)C)C3=C(OC)C=C(OC)C=C3C(=C2O)C4=C5OC(=CC(C5=C(O)C=6C4=CC(OC)=CC6OC)=O)C',
                                     'name': 'Isonigerone',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'C(=C\\C/C=C\\CCCC(NC)=O)\\C/C=C\\C/C=C\\CCCCC',
                                     'name': 'N-methyl arachidonoyl amine',
                                     'reason': 'No hydroxy groups found'},
                                 {   'smiles': 'OC1(O)[C@]23N(CC1)C(N(O)[C@H]([C@@]3(N=C(N2)N)[H])CO)=N',
                                     'name': 'Decarbamoylneosaxitoxin',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'S([C@H]1N(C(=O)/C(=C\\C2=CC=C(OCC=C(C)C)C=C2)/N(C1=O)C)C)C',
                                     'name': 'Fusaperazine F',
                                     'reason': 'No hydroxy groups found'},
                                 {   'smiles': 'Oc1ccccc1I',
                                     'name': '2-iodophenol',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'O(C1=C(OC)C=C(C2=C(O)C(OC)=C(C3=CC=CC=C3)C=C2OC)C=C1)C/C=C(/CO)\\C',
                                     'name': 'Prenylterphenyllin F',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)[C@@H](O)[C@@H]1OC[C@H]3O[C@@H](OC[C@H]4O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]4O)[C@H](O)[C@@H](O)[C@@H]3O)CO[C@@H]5O[C@@H]([C@@H](O)[C@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)[C@H]5O)CO[C@@H]7O[C@@H]([C@@H](O)[C@H](O)[C@H]7O)CO',
                                     'name': '(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5R,6R)-6-[[(2R,3R,4S,5R,6R)-3,5-Dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-2-yl]oxymethyl]-3,5-dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-3,4,5-trihydroxyoxan-2-yl]oxymethyl]oxane-2,3,4,5-tetrol',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'CN1CC2(CCN(CC2)S(=O)(=O)C3=CN(C=N3)C)C4=C([C@@H]1CO)NC5=C4C=CC(=C5)OC',
                                     'name': "[(1R)-7-methoxy-2-methyl-1'-[(1-methyl-4-imidazolyl)sulfonyl]-1-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,4'-piperidine]yl]methanol",
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'O=C(N(CCO)CCO)CN1C(=NC=C1)[N+]([O-])=O',
                                     'name': '1H-Imidazole-1-acetamide, '
                                             'N,N-bis(2-hydroxyethyl)-2-nitro-',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'CC1=NN(C(C1)(C2=CN=CC=C2)O)C(=O)CC3=CC=C(C=C3)C',
                                     'name': '1-[5-hydroxy-3-methyl-5-(3-pyridinyl)-4H-pyrazol-1-yl]-2-(4-methylphenyl)ethanone',
                                     'reason': 'No carboxylic acid group '
                                               'found'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To classify a chemical entity as a hydroxy fatty acid, we '
                 'need to verify two main features from the SMILES string: \n'
                 '\n'
                 '1. The presence of the carboxylic acid group `C(=O)O` or '
                 '`CO(=O)` which identifies the fatty acid backbone.\n'
                 '2. The existence of one or more hydroxyl groups `O`, '
                 'indicating the "hydroxy" characteristic of the molecule.\n'
                 '\n'
                 'Considering these two criteria, we can develop a function '
                 'using RDKit to parse the SMILES string and look for these '
                 'specific substructures in the given chemical. The function '
                 'will return `True` if both features are found and `False` '
                 'otherwise, along with a reason for the classification.\n'
                 '\n'
                 "Here's the program for this classification task:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 574,
    'num_false_positives': 100,
    'num_true_negatives': 90,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.8516320474777448,
    'recall': 1.0,
    'f1': 0.9198717948717948,
    'accuracy': 0.8691099476439791,
    'negative_predictive_value': 1.0}