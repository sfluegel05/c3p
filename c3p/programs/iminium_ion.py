"""
Classifies: CHEBI:35286 iminium ion
"""
from rdkit import Chem

def is_iminium_ion(smiles: str):
    """
    Determines if a molecule is an iminium ion (R2C=N(+)R2).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an iminium ion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for [N+]=C pattern
    pattern = Chem.MolFromSmarts('[N+]=C')
    if mol.HasSubstructMatch(pattern):
        return True, "Molecule contains the iminium ion structure R2C=N(+)R2"
    else:
        return False, "Molecule does not contain the iminium ion structure R2C=N(+)R2"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35286',
                          'name': 'iminium ion',
                          'definition': 'Cations of structure R2C=N(+)R2.',
                          'parents': ['CHEBI:25697']},
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
    'num_true_positives': 11,
    'num_false_positives': 0,
    'num_true_negatives': 11,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}