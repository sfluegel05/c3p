"""
Classifies: CHEBI:52575 alkylglycerol
"""
from rdkit import Chem

def is_alkylglycerol(smiles: str):
    """
    Determines if a molecule is an alkylglycerol (a glycerol ether having at least one alkyl substituent on oxygen).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkylglycerol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("[C@H](CO)(CO)O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for alkyl ether substituent on oxygen
    ether_pattern = Chem.MolFromSmarts("COC")
    if mol.HasSubstructMatch(ether_pattern):
        return True, "Molecule is an alkylglycerol"

    return False, "No alkyl ether substituent on oxygen found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:52575',
                          'name': 'alkylglycerol',
                          'definition': 'A glycerol ether having at least one '
                                        'alkyl substituent on oxygen.',
                          'parents': ['CHEBI:24353']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 48,
    'num_false_positives': 8,
    'num_true_negatives': 12,
    'num_false_negatives': 10,
    'precision': 0.8571428571428571,
    'recall': 0.8275862068965517,
    'f1': 0.8421052631578947,
    'accuracy': None}