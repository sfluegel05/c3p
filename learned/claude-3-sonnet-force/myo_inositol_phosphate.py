"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
"""
Classifies: CHEBI:27767 myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is an inositol phosphate in which the inositol component has myo-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for inositol backbone pattern
    inositol_pattern = Chem.MolFromSmarts("[C@@H]1([C@H]([C@@H]([C@H]([C@@H]([C@@H]1O)O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No myo-inositol backbone found"

    # Look for phosphate groups (-O-P(=O)(O)-O-)
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate groups found"

    # Count number of phosphate groups
    n_phosphates = len(phosphate_matches)

    # Optionally, you can add additional checks for chemical properties or molecular weight

    return True, f"Contains myo-inositol backbone with {n_phosphates} phosphate group(s)"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:27767',
        'name': 'myo-inositol phosphate',
        'definition': 'An inositol phosphate in which the inositol component has myo-configuration.',
        'parents': ['CHEBI:27761', 'CHEBI:26693']
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
    'num_true_positives': 170,
    'num_false_positives': 2,
    'num_true_negatives': 182409,
    'num_false_negatives': 3,
    'num_negatives': None,
    'precision': 0.9885057471264368,
    'recall': 0.9829545454545455,
    'f1': 0.9857142857142858,
    'accuracy': 0.9999619614771956
}