"""
Classifies: CHEBI:139592 tertiary alpha-hydroxy ketone
"""
from rdkit import Chem

def is_tertiary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a tertiary alpha-hydroxy ketone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for tertiary alpha-hydroxy ketone
    pattern = Chem.MolFromSmarts('[C](O)([C])[C](=O)[C]')
    if pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Molecule is a tertiary alpha-hydroxy ketone"
    else:
        return False, "Molecule does not match the tertiary alpha-hydroxy ketone pattern"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:139592',
                          'name': 'tertiary alpha-hydroxy ketone',
                          'definition': 'An alpha-hydroxy ketone in which the '
                                        'carbonyl group and the hydroxy group '
                                        'are linked by a carbon bearing two '
                                        'organyl groups.',
                          'parents': ['CHEBI:139588', 'CHEBI:26878']},
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
    'num_true_positives': 21,
    'num_false_positives': 16,
    'num_true_negatives': 4,
    'num_false_negatives': 1,
    'precision': 0.5675675675675675,
    'recall': 0.9545454545454546,
    'f1': 0.711864406779661,
    'accuracy': None}