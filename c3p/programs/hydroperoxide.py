"""
Classifies: CHEBI:35923 hydroperoxide
"""
from rdkit import Chem

def is_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a hydroperoxide (a monosubstitution product of hydrogen peroxide, HOOH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroperoxide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the hydroperoxide group pattern
    hydroperoxide_pattern = Chem.MolFromSmarts('OO')
    if hydroperoxide_pattern is None:
        return False, "Invalid hydroperoxide pattern"

    # Check if the molecule contains the hydroperoxide group
    if mol.HasSubstructMatch(hydroperoxide_pattern):
        return True, "Contains hydroperoxide group"
    else:
        return False, "Does not contain hydroperoxide group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35923',
                          'name': 'hydroperoxide',
                          'definition': 'A monosubstitution product of '
                                        'hydrogen peroxide, HOOH.',
                          'parents': ['CHEBI:24651']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 30-31: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}