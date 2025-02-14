"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
"""
Classifies: CHEBI:73199 2-hydroxydicarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid is a dicarboxylic acid carrying a hydroxy group on the carbon atom
    at position alpha to the carboxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxydicarboxylic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    carboxy_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    alpha_hydroxy_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H][CX4][OX2H]")

    # Find all carboxylic acid groups
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)
    num_carboxy_groups = len(carboxy_matches)
    if num_carboxy_groups < 2:
        return False, f"Found {num_carboxy_groups} carboxylic acid groups, need at least 2"

    # Find alpha-hydroxy carboxylic acid groups
    alpha_hydroxy_matches = mol.GetSubstructMatches(alpha_hydroxy_pattern)
    if not alpha_hydroxy_matches:
        return False, "No alpha-hydroxy carboxylic acid group found"

    return True, "Molecule is a 2-hydroxydicarboxylic acid (dicarboxylic acid with hydroxy group on alpha carbon)"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:73199',
        'name': '2-hydroxydicarboxylic acid',
        'definition': 'Any dicarboxylic acid carrying a hydroxy group on the carbon atom at position alpha to the carboxy group.',
        'parents': ['CHEBI:35692', 'CHEBI:14422']
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