"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: CHEBI:17410 quinic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is quinic acid or a derivative based on its SMILES string.
    Quinic acid is a cyclitol carboxylic acid with a cyclohexane backbone, multiple hydroxyl groups, and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is quinic acid or a derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the core cyclohexane ring with multiple hydroxyl groups and a carboxylic acid group
    quinic_acid_pattern = Chem.MolFromSmarts("[C@H]1([C@H](O)[C@@H](O)[C@H](O)[C@@H](C1)C(=O)O)")
    if not mol.HasSubstructMatch(quinic_acid_pattern):
        return False, "No quinic acid core structure found"

    # Count hydroxyl groups (OH) attached to the cyclohexane ring
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
    if hydroxyl_count < 3:
        return False, f"Found {hydroxyl_count} hydroxyl groups, need at least 3"

    # Check for the carboxylic acid group (C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for ester groups (common in quinic acid derivatives)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) > 0:
        return True, "Contains quinic acid core with ester groups (likely a derivative)"

    return True, "Contains quinic acid core structure with multiple hydroxyl groups and a carboxylic acid group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17410',
                          'name': 'quinic acid',
                          'definition': 'A cyclitol carboxylic acid.',
                          'parents': ['CHEBI:23403', 'CHEBI:23404']},
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