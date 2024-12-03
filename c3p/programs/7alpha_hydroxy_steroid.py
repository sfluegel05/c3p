"""
Classifies: CHEBI:36843 7alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_7alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 7alpha-hydroxy steroid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 7alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid backbone (cyclopentanoperhydrophenanthrene)
    steroid_backbone = Chem.MolFromSmarts('C1CCC2C1CCC3C2CCC4C3CCCC4')
    if not mol.HasSubstructMatch(steroid_backbone):
        return False, "No steroid backbone found"

    # Check for 7alpha-hydroxy group
    seven_alpha_hydroxy = Chem.MolFromSmarts('[C@H](O)[C@H](CCCC)C')
    if not mol.HasSubstructMatch(seven_alpha_hydroxy):
        return False, "No 7alpha-hydroxy group found"

    return True, "Molecule is a 7alpha-hydroxy steroid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36843',
                          'name': '7alpha-hydroxy steroid',
                          'definition': 'A 7-hydroxy steroid in which the '
                                        'hydroxy group at position 7 has an '
                                        'alpha-configuration.',
                          'parents': ['CHEBI:36844']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 10,
    'num_false_positives': 0,
    'num_true_negatives': 11,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9090909090909091,
    'f1': 0.9523809523809523,
    'accuracy': None}