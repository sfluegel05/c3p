"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: CHEBI:22680 arenecarbaldehyde
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is any aldehyde in which the carbonyl group is attached to an aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for aldehyde group (-CHO)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not aldehyde_matches:
        return False, "No aldehyde group found"

    # Check if the aldehyde carbon is directly attached to an aromatic atom
    aromatic_aldehyde_pattern = Chem.MolFromSmarts("[a]-[CX3H1]=O")
    aromatic_aldehyde_matches = mol.GetSubstructMatches(aromatic_aldehyde_pattern)
    if not aromatic_aldehyde_matches:
        return False, "Aldehyde group not directly attached to an aromatic atom"

    # Check if the molecule contains at least one aromatic ring
    aromatic_ring_pattern = Chem.MolFromSmarts("[a]")
    aromatic_ring_matches = mol.GetSubstructMatches(aromatic_ring_pattern)
    if not aromatic_ring_matches:
        return False, "No aromatic ring found in the molecule"

    # Ensure the molecule is primarily an aromatic aldehyde
    # by checking that the aldehyde is the only significant functional group
    # (excluding simple substituents like -OH, -OCH3, etc.)
    # This is a heuristic and may need refinement
    significant_functional_groups = ["[CX3](=O)[OX2H1]", "[NX3]", "[SX2]", "[PX3]"]
    for pattern in significant_functional_groups:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, f"Contains additional significant functional group: {pattern}"

    return True, "Contains an aldehyde group directly attached to an aromatic atom"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22680',
                          'name': 'arenecarbaldehyde',
                          'definition': 'Any aldehyde in which the carbonyl group is attached to an aromatic moiety.',
                          'parents': ['CHEBI:17478', 'CHEBI:22680']},
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