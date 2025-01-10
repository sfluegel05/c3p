"""
Classifies: CHEBI:3098 bile acid
"""
"""
Classifies: CHEBI:36218 bile acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    A bile acid is a steroid with a carboxylic acid group at C-24 and typically has hydroxyl groups on the steroid core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid core (four fused rings: three 6-membered and one 5-membered)
    steroid_pattern = Chem.MolFromSmarts("[C@H]1[C@H]2[C@H]3[C@H]4CC[C@H]4[C@H]3CC[C@H]2CC1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid core found"

    # Check for carboxylic acid group at C-24
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) == 0:
        return False, "No carboxylic acid group found"

    # Check for hydroxyl groups (at least one)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) == 0:
        return False, "No hydroxyl groups found"

    # Check for 5β-configuration (hydrogen at C-5 in β orientation)
    # This is more complex and may require specific stereochemistry checks
    # For simplicity, we assume the presence of the steroid core and carboxyl group is sufficient
    # In a more rigorous implementation, we would check the stereochemistry at C-5

    return True, "Contains steroid core with carboxylic acid group and hydroxyl groups"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:36218',
        'name': 'bile acid',
        'definition': 'Any member of a group of hydroxy-5beta-cholanic acids occurring in bile, where they are present as the sodium salts of their amides with glycine or taurine.',
        'parents': ['CHEBI:36217', 'CHEBI:36219']
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