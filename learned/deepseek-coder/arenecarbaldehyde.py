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

    # Check if the aldehyde carbon is attached to an aromatic atom
    aromatic_aldehyde_pattern = Chem.MolFromSmarts("[c;$(c-[CX3H1]=O)]")
    aromatic_aldehyde_matches = mol.GetSubstructMatches(aromatic_aldehyde_pattern)
    if not aromatic_aldehyde_matches:
        return False, "Aldehyde group not attached to an aromatic ring"

    # Check if the molecule contains at least one aromatic ring
    aromatic_ring_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1")
    aromatic_ring_matches = mol.GetSubstructMatches(aromatic_ring_pattern)
    if not aromatic_ring_matches:
        return False, "No aromatic ring found in the molecule"

    return True, "Contains an aldehyde group attached to an aromatic ring"


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