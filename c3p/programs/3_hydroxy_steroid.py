"""
Classifies: CHEBI:36834 3-hydroxy steroid
"""
from rdkit import Chem

def is_3_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy steroid (any hydroxy steroid carrying a hydroxy group at position 3).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a steroid backbone
    steroid_backbone = Chem.MolFromSmarts("C1CCC2C1CCC3C2CCC4C3CCC(C4)C")
    if not mol.HasSubstructMatch(steroid_backbone):
        return False, "No steroid backbone found"

    # Check for the presence of a hydroxy group at position 3
    hydroxy_at_3 = Chem.MolFromSmarts("C1[C@H](O)CC2C1CCC3C2CCC4C3CCC(C4)C")
    if mol.HasSubstructMatch(hydroxy_at_3):
        return True, "3-hydroxy steroid found"
    
    # Check for the presence of a hydroxy group at position 3 (alternative stereochemistry)
    hydroxy_at_3_alt = Chem.MolFromSmarts("C1[C@@H](O)CC2C1CCC3C2CCC4C3CCC(C4)C")
    if mol.HasSubstructMatch(hydroxy_at_3_alt):
        return True, "3-hydroxy steroid found"

    return False, "No hydroxy group at position 3 found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36834',
                          'name': '3-hydroxy steroid',
                          'definition': 'Any hydroxy steroid carrying a '
                                        'hydroxy group at position 3.',
                          'parents': ['CHEBI:35350']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[23:58:40] SMILES Parse Error: syntax error while parsing: '
             'O=C\x01C[C@]2([C@]3([C@@]([C@@]4(C(CC3)=CC(=O)CC4)C)(CC[C@@]2(/C1=C\\C)C)[H])[H])[H]\n'
             '[23:58:40] SMILES Parse Error: Failed parsing SMILES '
             "'O=C\x01C[C@]2([C@]3([C@@]([C@@]4(C(CC3)=CC(=O)CC4)C)(CC[C@@]2(/C1=C\\C)C)[H])[H])[H]' "
             'for input: '
             "'O=C\x01C[C@]2([C@]3([C@@]([C@@]4(C(CC3)=CC(=O)CC4)C)(CC[C@@]2(/C1=C\\C)C)[H])[H])[H]'\n",
    'stdout': '',
    'num_true_positives': 1,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 83,
    'precision': 1.0,
    'recall': 0.011904761904761904,
    'f1': 0.023529411764705882,
    'accuracy': None}