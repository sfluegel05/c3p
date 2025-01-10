"""
Classifies: CHEBI:22307 aldoxime
"""
"""
Classifies: CHEBI:36586 aldoxime
"""
from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime is characterized by the presence of the RCH=NOH group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldoxime, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the aldoxime pattern: R-CH=NOH
    aldoxime_pattern = Chem.MolFromSmarts("[CX3H1](=N[OH])")
    if not mol.HasSubstructMatch(aldoxime_pattern):
        return False, "No aldoxime group (R-CH=NOH) found"

    # Check if the carbon in the aldoxime group is part of an aldehyde
    # An aldehyde carbon is connected to a hydrogen and another carbon or heteroatom
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No aldehyde group (R-CHO) found"

    return True, "Contains the aldoxime group (R-CH=NOH) derived from an aldehyde"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36586',
                          'name': 'aldoxime',
                          'definition': 'Oximes of aldehydes RCH=NOH.',
                          'parents': ['CHEBI:36585', 'CHEBI:36587']},
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