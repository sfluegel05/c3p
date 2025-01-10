"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: essential fatty acid
Any member of the sub-set of polyunsaturated fatty acid for which there is an absolute dietary requirement.
"""
def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: None
        str: None
    """
    return None, None

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'essential fatty acid',
        'definition': 'Any member of the sub-set of polyunsaturated fatty acid for which there is an absolute dietary requirement.',
        'parents': ['polyunsaturated fatty acid', 'fatty acid'],
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
    'attempt': 1,
    'success': False,
    'best': False,
    'error': '',
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}