"""
Classifies: CHEBI:36844 7-hydroxy steroid
"""
from rdkit import Chem

def is_7_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 7-hydroxy steroid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 7-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains a steroid backbone (cyclopentanoperhydrophenanthrene)
    steroid_scaffold = Chem.MolFromSmarts('C1CCC2C1CCC3C2CCC4C3CCCC4')
    if not mol.HasSubstructMatch(steroid_scaffold):
        return False, "Molecule does not contain steroid backbone"

    # Check for 7-hydroxy substituent
    hydroxy_7 = Chem.MolFromSmarts('[C@@H]1CC[C@@]2(C)[C@]3([H])CC[C@]4(C)[C@]([H])(CC[C@@]4([H])[C@]3([H])[C@H](O)C2)[C@H](C)C1')
    if mol.HasSubstructMatch(hydroxy_7):
        return True, "Molecule is a 7-hydroxy steroid"
    else:
        return False, "Molecule does not have a hydroxy substituent at position 7"

# Example usage:
# smiles = "C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3([C@@H](C[C@@]2(C[C@H](C1)O)[H])O)[H])(CC[C@@]4([C@@H](CCC(O)=O)C)[H])[H])C)[H])C"
# result, reason = is_7_hydroxy_steroid(smiles)
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36844',
                          'name': '7-hydroxy steroid',
                          'definition': 'A hydroxy steroid carrying a hydroxy '
                                        'substituent at position 7.',
                          'parents': ['CHEBI:35350']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[00:00:59] SMILES Parse Error: syntax error while parsing: '
             'Cl/C=C\x01/C=CC=2C=C3C=C(CC(=O)O)OC3=CC2OC1\n'
             '[00:00:59] SMILES Parse Error: Failed parsing SMILES '
             "'Cl/C=C\x01/C=CC=2C=C3C=C(CC(=O)O)OC3=CC2OC1' for input: "
             "'Cl/C=C\x01/C=CC=2C=C3C=C(CC(=O)O)OC3=CC2OC1'\n",
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 14,
    'num_false_negatives': 14,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}