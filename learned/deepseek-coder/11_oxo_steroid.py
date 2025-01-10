"""
Classifies: CHEBI:47787 11-oxo steroid
"""
"""
Classifies: CHEBI:78298 11-oxo steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    An 11-oxo steroid is a steroid with an oxo group (C=O) at position 11.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11-oxo steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more specific steroid backbone pattern (four fused rings)
    steroid_backbone = Chem.MolFromSmarts("[C@]12[C@H]3[C@H]([C@H]4[C@H]([C@H]1CC2)CC3)CC4")
    if not mol.HasSubstructMatch(steroid_backbone):
        return False, "No steroid backbone found"

    # Define the oxo group at position 11
    oxo_at_11_pattern = Chem.MolFromSmarts("[C@@H]1[C@H]2[C@H]([C@H]3[C@H]([C@H]1CC2)CC3)C(=O)")
    if mol.HasSubstructMatch(oxo_at_11_pattern):
        return True, "Contains steroid backbone with oxo group at position 11"

    return False, "No oxo group found at position 11"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:78298',
                          'name': '11-oxo steroid',
                          'definition': 'Any oxo steroid that has an oxo substituent at position 11.',
                          'parents': ['CHEBI:78298', 'CHEBI:78298']},
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
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}