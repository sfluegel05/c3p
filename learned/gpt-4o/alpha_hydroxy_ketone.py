"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
from rdkit import Chem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is a ketone containing a hydroxy group on the alpha-carbon 
    relative to the C=O group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ketone pattern: O=C-C(alpha)
    ketone_pattern = Chem.MolFromSmarts("C(=O)C")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group found"

    # Alpha-hydroxy pattern: O=C-C(O)
    alpha_hydroxy_pattern = Chem.MolFromSmarts("C(=O)C(O)")
    if mol.HasSubstructMatch(alpha_hydroxy_pattern):
        return True, "Contains alpha-hydroxy ketone pattern"
    else:
        return False, "No hydroxy group on alpha-carbon"

__metadata__ = {   'chemical_class': {   'id': None,
                                         'name': 'alpha-hydroxy ketone',
                                         'definition': 'A ketone containing a hydroxyl group on the alpha-carbon relative to the C=O group.'},
                    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                                  'f1_threshold': 0.8,
                                  'max_attempts': 5,
                                  'max_positive_instances': None,
                                  'max_positive_to_test': None,
                                  'max_negative_to_test': None,
                                  'max_positive_in_prompt': 50,
                                  'max_negative_in_prompt': 20,
                                  'max_instances_in_prompt': 100,
                                  'test_proportion': 0.1},
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
                    'precision': None,
                    'recall': None,
                    'f1': None,
                    'accuracy': None}