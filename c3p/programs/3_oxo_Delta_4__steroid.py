"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
from rdkit import Chem


def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(4) steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a ketone group at position 3 (3-oxo)
    ketone_3 = Chem.MolFromSmarts('C(=O)')
    matches = mol.GetSubstructMatches(ketone_3)
    if not any(match[0] == 2 for match in matches):  # Position 3 in a steroid scaffold
        return False, "No 3-oxo group found"

    # Check for a double bond at the Delta(4) position (C=C at position 4)
    double_bond_4 = Chem.MolFromSmarts('C=CC')
    matches = mol.GetSubstructMatches(double_bond_4)
    if not any(match[0] == 3 for match in matches):  # Position 4 in a steroid scaffold
        return False, "No Delta(4) double bond found"

    # Check for the steroid scaffold
    steroid_scaffold = Chem.MolFromSmarts('C1CCC2C1CCC3C2CCC4C3CC=C4')
    if not mol.HasSubstructMatch(steroid_scaffold):
        return False, "No steroid scaffold found"

    return True, "Molecule is a 3-oxo-Delta(4) steroid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:47909',
                          'name': '3-oxo-Delta(4) steroid',
                          'definition': 'A 3-oxo steroid conjugated to a C=C '
                                        'double bond at the alpha,beta '
                                        'position.',
                          'parents': ['CHEBI:47788', 'CHEBI:51689']},
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
    'num_true_negatives': 20,
    'num_false_negatives': 25,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}