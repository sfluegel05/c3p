"""
Classifies: CHEBI:74716 withanolide
"""
from rdkit import Chem

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid backbone (C28 steroid)
    steroid_pattern = Chem.MolFromSmarts('[C]1([C@H]2[C@H]3[C@H]4[C@H]5[C@@H]6[C@H]7[C@@H]8[C@@H]9[C@@H]1[C@@H]2[C@@H]3[C@@H]4[C@@H]5[C@@H]6[C@@H]7[C@@H]8[C@@H]9)')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for lactone ring
    lactone_pattern = Chem.MolFromSmarts('O=C1OC[C@@H]1')
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"

    # Check for modified side chain forming a lactone ring
    side_chain_lactone_pattern = Chem.MolFromSmarts('O=C1O[C@@H]([C@H]2[C@@H]([C@@H]3[C@@H]([C@H]4[C@@H]([C@@H](C1)C)C)C)C)C')
    if not mol.HasSubstructMatch(side_chain_lactone_pattern):
        return False, "No modified side chain forming a lactone ring found"

    return True, "Molecule is a withanolide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:74716',
                          'name': 'withanolide',
                          'definition': 'Any steroid lactone that is a C28 '
                                        'steroid with a modified side chain '
                                        'forming a lactone ring and its '
                                        'substituted derivatives.',
                          'parents': ['CHEBI:26766']},
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
    'num_true_negatives': 15,
    'num_false_negatives': 15,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}