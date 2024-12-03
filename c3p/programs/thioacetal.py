"""
Classifies: CHEBI:59792 thioacetal
"""
from rdkit import Chem

def is_thioacetal(smiles: str):
    """
    Determines if a molecule is a thioacetal.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thioacetal, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for monothioacetal and dithioacetal
    monothioacetal_smarts = Chem.MolFromSmarts("[#6]([#8][#6])([#16][#6])")
    dithioacetal_smarts = Chem.MolFromSmarts("[#6]([#16][#6])([#16][#6])")

    if mol.HasSubstructMatch(monothioacetal_smarts):
        return True, "Monothioacetal structure found"
    elif mol.HasSubstructMatch(dithioacetal_smarts):
        return True, "Dithioacetal structure found"
    else:
        return False, "No thioacetal structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:59792',
                          'name': 'thioacetal',
                          'definition': "The sulfur analogue of 'acetal'.  The "
                                        'term includes monothioacetals having '
                                        "the structure R2C(OR')(SR') (subclass "
                                        'monothioketals, R =/= H); and '
                                        'dithioacetals having the structure '
                                        "R2C(SR')2 (subclass dithioketals, R "
                                        "=/= H, R' =/= H).",
                          'parents': ['CHEBI:33261']},
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
    'num_true_negatives': 10,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}