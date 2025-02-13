"""
Classifies: CHEBI:15734 primary alcohol
"""
from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the primary alcohol pattern: -CH2OH
    primary_alcohol_pattern = Chem.MolFromSmarts("[CH2]O")

    # Check for primary alcohol pattern
    if mol.HasSubstructMatch(primary_alcohol_pattern):
        return True, "Contains primary alcohol group"

    return False, "Does not contain primary alcohol group"


__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'primary alcohol',
        'definition': 'A primary alcohol is a compound in which a hydroxy group, -OH, is attached to a saturated carbon atom which has either three hydrogen atoms attached to it or only one other carbon atom and two hydrogen atoms attached to it.',
        'parents': []
    },
    'config': {   
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0,
    'recall': 0,
    'f1': 0,
    'accuracy': 0
}