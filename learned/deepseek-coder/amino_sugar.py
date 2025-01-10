"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: CHEBI:37671 amino sugar
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is a sugar with one or more hydroxyl groups replaced by amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is a sugar (contains multiple hydroxyl groups and a ring)
    sugar_pattern = Chem.MolFromSmarts("[OH]-[C]-[C]-[C]-[C]-[C]-[OH]")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar-like structure found"

    # Check for the presence of amino groups (NH2 or substituted amino groups)
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    if len(amino_matches) == 0:
        return False, "No amino groups found"

    # Check if the amino group is attached to a carbon in the sugar ring
    amino_sugar_pattern = Chem.MolFromSmarts("[C]-[NX3;H2,H1;!$(NC=O)]")
    if not mol.HasSubstructMatch(amino_sugar_pattern):
        return False, "Amino group not attached to sugar carbon"

    # Count the number of amino groups
    n_amino = len(amino_matches)
    if n_amino < 1:
        return False, "No amino groups found"

    # Count the number of hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    n_hydroxyl = len(hydroxyl_matches)

    # Ensure the molecule has multiple hydroxyl groups (sugar-like)
    if n_hydroxyl < 2:
        return False, "Not enough hydroxyl groups for a sugar"

    return True, f"Contains a sugar structure with {n_amino} amino group(s) attached"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37671',
                          'name': 'amino sugar',
                          'definition': 'Any sugar having one or more alcoholic '
                                        'hydroxy groups replaced by substituted '
                                        'or unsubstituted amino groups.',
                          'parents': ['CHEBI:47778', 'CHEBI:76579']},
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