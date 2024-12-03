"""
Classifies: CHEBI:24632 hydrocarbon
"""
from rdkit import Chem

def is_hydrocarbon(smiles: str):
    """
    Determines if a molecule is a hydrocarbon (a compound consisting of carbon and hydrogen only).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydrocarbon, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ['C', 'H']:
            return False, f"Contains non-carbon and non-hydrogen atom: {atom.GetSymbol()}"

    return True, "Molecule is a hydrocarbon"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24632',
                          'name': 'hydrocarbon',
                          'definition': 'A compound consisting of carbon and '
                                        'hydrogen only.',
                          'parents': ['CHEBI:33245']},
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
    'num_true_positives': 107,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 43,
    'precision': 1.0,
    'recall': 0.7133333333333334,
    'f1': 0.8326848249027238,
    'accuracy': None}