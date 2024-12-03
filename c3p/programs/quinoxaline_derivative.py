"""
Classifies: CHEBI:38771 quinoxaline derivative
"""
from rdkit import Chem

def is_quinoxaline_derivative(smiles: str):
    """
    Determines if a molecule is a quinoxaline derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinoxaline derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the quinoxaline core (1,4-naphthyridine)
    quinoxaline_core = Chem.MolFromSmarts('c1cnc2nccnc2c1')

    if quinoxaline_core is None:
        return False, "Error in defining quinoxaline core"

    # Check if the molecule contains the quinoxaline core
    if not mol.HasSubstructMatch(quinoxaline_core):
        return False, "Molecule does not contain the quinoxaline core"

    return True, "Molecule is a quinoxaline derivative"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38771',
                          'name': 'quinoxaline derivative',
                          'definition': 'Any naphthyridine derivative that is '
                                        'a derivative of quinoxaline '
                                        '(1,4-naphthyridine).',
                          'parents': ['CHEBI:73539']},
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 15,
    'num_false_negatives': 15,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}