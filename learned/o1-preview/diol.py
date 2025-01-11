"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol is a compound that contains two hydroxy groups, generally assumed to be, but not necessarily, alcoholic.
    Aliphatic diols are also called glycols.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for hydroxyl groups (OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    hydroxyl_count = len(hydroxyl_matches)

    if hydroxyl_count >= 2:
        return True, f"Molecule contains at least two hydroxyl groups"
    else:
        return False, f"Molecule contains {hydroxyl_count} hydroxyl group(s), diols must have at least two"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23824',
                              'name': 'diol',
                              'definition': 'A compound that contains two hydroxy '
                                            'groups, generally assumed to be, but not '
                                            'necessarily, alcoholic. Aliphatic diols are '
                                            'also called glycols.',
                              'parents': []},
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
        'attempt': 2,
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
        'accuracy': None}