"""
Classifies: CHEBI:59768 aliphatic aldehyde
"""
from rdkit import Chem

def is_aliphatic_aldehyde(smiles: str):
    """
    Determines if a molecule is an aliphatic aldehyde.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic aldehyde, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains an aldehyde group
    aldehyde = Chem.MolFromSmarts('[CX3H1](=O)[#6]')
    if not mol.HasSubstructMatch(aldehyde):
        return False, "No aldehyde group found"

    # Check if the molecule is aliphatic (contains no aromatic rings)
    if mol.GetRingInfo().NumAromaticRings() > 0:
        return False, "Contains aromatic rings"

    return True, "Aliphatic aldehyde"

# Example usage:
# smiles = "O=CCCCCCCC/C=C\CCCCCC"  # 9Z-Hexadecenal
# result, reason = is_aliphatic_aldehyde(smiles)
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:59768',
                          'name': 'aliphatic aldehyde',
                          'definition': 'An aldehyde derived from an aliphatic '
                                        'compound',
                          'parents': ['CHEBI:17478']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "'RingInfo' object has no attribute 'NumAromaticRings'",
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}