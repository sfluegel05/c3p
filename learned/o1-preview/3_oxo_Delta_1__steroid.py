"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
"""
Classifies: CHEBI:35164 3-oxo-Delta(1) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid based on its SMILES string.
    A 3-oxo-Delta(1) steroid is any 3-oxo steroid that contains a double bond between positions 1 and 2.
    This function checks for:
      - Steroid backbone (6-6-6-5 fused ring system)
      - Ketone group at position 3
      - Double bond between positions 1 and 2

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(1) steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid backbone pattern (6-6-6-5 fused ring system)
    steroid_pattern = Chem.MolFromSmarts("""
    [#6]1CC[C@H]2[C@H](C1)CC[C@@H]3[C@@H](C2)CC[C@@H]4C3=CC=CC4
    """)
    if steroid_pattern is None:
        return False, "Invalid steroid SMARTS pattern"

    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for ketone group at position 3
    ketone_pattern = Chem.MolFromSmarts("""
    [#6]1=CC[C@H]2[C@H](C1)CC[C@@H]3[C@@H](C2)CC[C@@H]4C3=CC(=O)C=C4
    """)
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group at position 3"

    # Check for double bond between positions 1 and 2
    double_bond_pattern = Chem.MolFromSmarts("""
    [#6]=1C=CC[C@H]2[C@H](C1)CC[C@@H]3[C@@H](C2)CC[C@@H]4C3=CC(=O)C=C4
    """)
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No double bond between positions 1 and 2"

    return True, "Matches 3-oxo-Delta(1) steroid with correct backbone, ketone at position 3, and double bond between positions 1 and 2"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35164',
                              'name': '3-oxo-Delta(1) steroid',
                              'definition': 'Any 3-oxo steroid that contains a '
                                            'double bond between positions 1 and 2.'},
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
        'stdout': None}