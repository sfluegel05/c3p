"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: CHEBI:32955 epoxide
"""
from rdkit import Chem

def is_epoxide(smiles: str):
    """
    Determines if a molecule is an epoxide based on its SMILES string.
    An epoxide is a cyclic ether in which the oxygen atom forms part of a 3-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the epoxide (oxirane) SMARTS pattern
    epoxide_pattern = Chem.MolFromSmarts("C1OC1")  # 3-membered ring with two carbons and one oxygen
    if epoxide_pattern is None:
        return False, "Invalid SMARTS pattern for epoxide"
    
    # Check for epoxide substructure
    if mol.HasSubstructMatch(epoxide_pattern):
        return True, "Contains an epoxide functional group (3-membered cyclic ether)"
    else:
        return False, "No epoxide functional group found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:32955',
        'name': 'epoxide',
        'definition': 'Any cyclic ether in which the oxygen atom forms part of a 3-membered ring.',
        'parents': ['CHEBI:23816', 'CHEBI:62866']
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
    # Placeholder metrics - these would be populated during testing
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