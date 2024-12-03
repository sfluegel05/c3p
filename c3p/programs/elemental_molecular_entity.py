"""
Classifies: CHEBI:33259 elemental molecular entity
"""
from rdkit import Chem

def is_elemental_molecular_entity(smiles: str):
    """
    Determines if a molecule is an elemental molecular entity (all atoms have the same atomic number).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an elemental molecular entity, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    atoms = mol.GetAtoms()
    if not atoms:
        return False, "No atoms found in the molecule"

    atomic_numbers = {atom.GetAtomicNum() for atom in atoms}

    if len(atomic_numbers) == 1:
        return True, f"All atoms have the same atomic number: {atomic_numbers.pop()}"
    else:
        return False, f"Atoms have different atomic numbers: {', '.join(map(str, atomic_numbers))}"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33259',
                          'name': 'elemental molecular entity',
                          'definition': 'A molecular entity all atoms of which '
                                        'have the same atomic number.',
                          'parents': ['CHEBI:23367']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[22:59:34] Explicit valence for atom # 1 H, 2, is greater than '
             'permitted\n'
             '[22:59:34] Explicit valence for atom # 1 Br, 2, is greater than '
             'permitted\n'
             '[22:59:34] Explicit valence for atom # 0 B, 5, is greater than '
             'permitted\n'
             '[22:59:34] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[22:59:34] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[22:59:34] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[22:59:34] SMILES Parse Error: syntax error while parsing: '
             'O(C([C@H](\\C=C\\[C@H]([C@@]1([C@@]2(C(CC1)/C(/CCC2)=C/C=C\x03/CC(O)CCC3=C)C)[H])C)C)(C)C)[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)C(O)=O\n'
             '[22:59:34] SMILES Parse Error: Failed parsing SMILES '
             "'O(C([C@H](\\C=C\\[C@H]([C@@]1([C@@]2(C(CC1)/C(/CCC2)=C/C=C\x03/CC(O)CCC3=C)C)[H])C)C)(C)C)[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)C(O)=O' "
             'for input: '
             "'O(C([C@H](\\C=C\\[C@H]([C@@]1([C@@]2(C(CC1)/C(/CCC2)=C/C=C\x03/CC(O)CCC3=C)C)[H])C)C)(C)C)[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)C(O)=O'\n",
    'stdout': '',
    'num_true_positives': 25,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 3,
    'precision': 0.9615384615384616,
    'recall': 0.8928571428571429,
    'f1': 0.9259259259259259,
    'accuracy': None}