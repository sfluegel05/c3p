"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: CHEBI:132849 arenecarbaldehyde
"""
from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is any aldehyde in which the carbonyl group is attached to an aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define aldehyde group attached to aromatic ring
    arenecarbaldehyde_pattern = Chem.MolFromSmarts("[a][CX3H1](=O)")

    # Check for aldehyde group attached to aromatic ring
    if not mol.HasSubstructMatch(arenecarbaldehyde_pattern):
        return False, "No aldehyde group attached to aromatic ring found"

    return True, "Contains aldehyde group attached to aromatic moiety"


__metadata__ = {   
    'chemical_class': {   
        'id': 'CHEBI:132849',
        'name': 'arenecarbaldehyde',
        'definition': 'Any aldehyde in which the carbonyl group is attached to an aromatic moiety.',
        'parents': ['CHEBI:17480', 'CHEBI:33690']
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
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199
}