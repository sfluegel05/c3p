"""
Classifies: CHEBI:64482 phosphatidylcholine
"""
"""
Classifies: phosphatidylcholine
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylcholine(smiles: str):
    """
    Determines if a molecule is a phosphatidylcholine based on its SMILES string.
    A phosphatidylcholine is a glycerophosphocholine that is glycero-3-phosphocholine bearing
    two acyl substituents at positions 1 and 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylcholine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphocholine head group
    phosphocholine_smarts = "[P](=O)([O-])(OCC[N+](C)(C)C)"
    phosphocholine_pattern = Chem.MolFromSmarts(phosphocholine_smarts)
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine head group found"

    # Check for glycerol backbone with two ester groups
    glycerol_esters_smarts = "[C@@H](CO[P](=O)([O-])OCC[N+](C)(C)C)(OC(=O)[#6])OC(=O)[#6]"
    glycerol_esters_pattern = Chem.MolFromSmarts(glycerol_esters_smarts)
    if not mol.HasSubstructMatch(glycerol_esters_pattern):
        return False, "No glycerol backbone with two ester groups found"

    return True, "Contains glycerol backbone with two acyl chains and phosphocholine head group"

__metadata__ = {
   'chemical_class': {
      'name': 'phosphatidylcholine',
      'definition': 'A glycerophosphocholine that is glycero-3-phosphocholine bearing two acyl substituents at positions 1 and 2.',
      'id': None,
      'parents': []
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