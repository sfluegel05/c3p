"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for 17beta-hydroxy group
    pattern = Chem.MolFromSmarts('[C@H](O)[C@@]([C@H]1CC[C@@H]2[C@@H](CCC3=CC(=O)CC[C@]32C)[C@@H]1C)C')
    if pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check if the molecule matches the 17beta-hydroxy steroid pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Molecule is a 17beta-hydroxy steroid"
    else:
        return False, "Molecule does not match 17beta-hydroxy steroid pattern"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35343',
                          'name': '17beta-hydroxy steroid',
                          'definition': 'A 17-hydroxy steroid in which the '
                                        'hydroxy group at position 17 has a '
                                        'beta-configuration.',
                          'parents': ['CHEBI:36838']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 11,
    'num_false_negatives': 11,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}