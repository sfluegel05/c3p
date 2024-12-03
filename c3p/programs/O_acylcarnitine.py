"""
Classifies: CHEBI:17387 O-acylcarnitine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_O_acylcarnitine(smiles: str):
    """
    Determines if a molecule is an O-acylcarnitine (Any carboxylic ester obtained by the O-acylation of carnitine).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an O-acylcarnitine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the carnitine backbone: [C@H](CC([O-])=O)C[N+](C)(C)C
    carnitine_pattern = Chem.MolFromSmarts('[C@H](CC([O-])=O)C[N+](C)(C)C')
    if not mol.HasSubstructMatch(carnitine_pattern):
        return False, "Carnitine backbone not found"

    # Check for the ester linkage: O=C-O
    ester_pattern = Chem.MolFromSmarts('O=C-O')
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "Ester linkage not found"

    # Check that the ester linkage is connected to the carnitine backbone
    carnitine_matches = mol.GetSubstructMatches(carnitine_pattern)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    for c_match in carnitine_matches:
        for e_match in ester_matches:
            # Check if the ester oxygen is connected to the carnitine carbon
            if mol.GetBondBetweenAtoms(c_match[0], e_match[2]) is not None:
                return True, "O-acylcarnitine structure found"

    return False, "Ester linkage not connected to carnitine backbone"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17387',
                          'name': 'O-acylcarnitine',
                          'definition': 'Any carboxylic ester obtained by the '
                                        'O-acylation of carnitine.',
                          'parents': ['CHEBI:33308', 'CHEBI:35284']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 27,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 7,
    'precision': 1.0,
    'recall': 0.7941176470588235,
    'f1': 0.8852459016393442,
    'accuracy': None}