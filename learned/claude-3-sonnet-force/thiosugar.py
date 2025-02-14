"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: CHEBI:35656 thiosugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar is a carbohydrate derivative in which one or more of the
    oxygens or hydroxy groups of the parent carbohydrate is replaced by
    sulfur or -SR, where R can be hydrogen or any group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Detect carbohydrate backbone
    carbohydrate_pattern = Chem.MolFromSmarts("[OX2]([CX4]([OX2])([CX4]([OX2])([CX4]([OX2])))([H])[H]")
    if not mol.HasSubstructMatch(carbohydrate_pattern):
        return False, "No carbohydrate backbone found"

    # Detect sulfur replacing oxygen or hydroxy group
    thiosugar_pattern = Chem.MolFromSmarts("[SX2]([CX4]([OX2])([CX4]([OX2])([CX4]([OX2])))([H])[H]")
    if not mol.HasSubstructMatch(thiosugar_pattern):
        return False, "No sulfur atom replacing oxygen or hydroxy group"

    return True, "Contains a carbohydrate backbone with sulfur replacing oxygen or hydroxy group"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35656',
        'name': 'thiosugar',
        'definition': 'A carbohydrate derivative in which one or more of the oxygens or hydroxy groups of the parent carbohydrate is replaced by sulfur or -SR, where R can be hydrogen or any group.',
        'parents': ['CHEBI:24313', 'CHEBI:33567']
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
    'num_true_positives': 33,
    'num_false_positives': 0,
    'num_true_negatives': 182382,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': 0.9998192006472493
}