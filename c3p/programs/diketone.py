"""
Classifies: CHEBI:46640 diketone
"""
from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone (contains two ketone functionalities).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diketone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the ketone functional group pattern
    ketone_pattern = Chem.MolFromSmarts('C(=O)')

    # Find all ketone matches in the molecule
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    # Check if there are at least two ketone functionalities
    if len(ketone_matches) >= 2:
        return True, "Contains two or more ketone functionalities"
    else:
        return False, "Less than two ketone functionalities found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:46640',
                          'name': 'diketone',
                          'definition': 'A compound that contains two ketone '
                                        'functionalities.',
                          'parents': ['CHEBI:17087']},
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
    'num_true_positives': 18,
    'num_false_positives': 8,
    'num_true_negatives': 11,
    'num_false_negatives': 1,
    'precision': 0.6923076923076923,
    'recall': 0.9473684210526315,
    'f1': 0.7999999999999999,
    'accuracy': None}