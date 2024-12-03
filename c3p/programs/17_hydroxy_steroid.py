"""
Classifies: CHEBI:36838 17-hydroxy steroid
"""
from rdkit import Chem

def is_17_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17-hydroxy steroid (a hydroxy steroid carrying a hydroxy group at position 17).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is a steroid
    steroid_core_smarts = '[#6]1[#6][#6][#6]2[#6][#6][#6]3[#6][#6][#6]4[#6][#6][#6]([#6][#6]([#6]4[#6]3[#6]2[#6]1)[#6])'
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    if not mol.HasSubstructMatch(steroid_core):
        return False, "Molecule does not contain the steroid core structure"

    # Check for hydroxy group at position 17
    hydroxy_17_smarts = '[#6]1[#6][#6][#6]2[#6][#6][#6]3[#6][#6][#6]4[#6][#6][#6]([#6][#6]([#6]4[#6]3[#6]2[#6]1)[#6])[#8H]'
    hydroxy_17 = Chem.MolFromSmarts(hydroxy_17_smarts)
    if not mol.HasSubstructMatch(hydroxy_17):
        return False, "Molecule does not have a hydroxy group at position 17"

    return True, "Molecule is a 17-hydroxy steroid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36838',
                          'name': '17-hydroxy steroid',
                          'definition': 'A hydroxy steroid carrying a hydroxy '
                                        'group at position 17.',
                          'parents': ['CHEBI:35350']},
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
    'num_false_negatives': 24,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}