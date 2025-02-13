"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: CHEBI:50584 aliphatic alcohol
"""
from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol is defined as an alcohol derived from an aliphatic compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one hydroxyl group (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # Check that the hydroxyl group is attached to an aliphatic carbon
    aliphatic_hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OH]")
    if not mol.HasSubstructMatch(aliphatic_hydroxyl_pattern):
        return False, "Hydroxyl group not attached to an aliphatic carbon"

    # Ensure the molecule is not aromatic (no aromatic rings)
    if mol.GetAromaticAtoms():
        return False, "Molecule contains aromatic rings"

    return True, "Contains a hydroxyl group attached to an aliphatic carbon chain"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50584',
                          'name': 'aliphatic alcohol',
                          'definition': 'An alcohol derived from an aliphatic compound.',
                          'parents': ['CHEBI:23888', 'CHEBI:50584']},
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