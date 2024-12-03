"""
Classifies: CHEBI:26588 1,3,5-triazines
"""
from rdkit import Chem

def is_1_3_5_triazines(smiles: str):
    """
    Determines if a molecule is a 1,3,5-triazine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1,3,5-triazine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for 1,3,5-triazine
    triazine_pattern = Chem.MolFromSmarts("n1cncnc1")

    if mol.HasSubstructMatch(triazine_pattern):
        return True, "Molecule contains the 1,3,5-triazine skeleton"
    else:
        return False, "Molecule does not contain the 1,3,5-triazine skeleton"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26588',
                          'name': '1,3,5-triazines',
                          'definition': 'Any compound with a 1,3,5-triazine '
                                        'skeleton, in which nitrogen atoms '
                                        'replace carbon at positions 1, 3 and '
                                        '5 of the core benzene ring structure.',
                          'parents': ['CHEBI:38102']},
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
    'num_true_positives': 13,
    'num_false_positives': 1,
    'num_true_negatives': 12,
    'num_false_negatives': 0,
    'precision': 0.9285714285714286,
    'recall': 1.0,
    'f1': 0.962962962962963,
    'accuracy': None}