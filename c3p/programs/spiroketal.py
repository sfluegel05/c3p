"""
Classifies: CHEBI:72600 spiroketal
"""
from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a spiroketal, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    rings = mol.GetRingInfo().AtomRings()

    # Identify spiro atoms (atoms shared by two rings)
    spiro_atoms = set()
    for i, ring1 in enumerate(rings):
        for ring2 in rings[i+1:]:
            common_atoms = set(ring1) & set(ring2)
            if len(common_atoms) == 1:
                spiro_atoms.update(common_atoms)

    if not spiro_atoms:
        return False, "No spiro atoms found"

    # Check if the spiro atoms are ketal carbons
    for atom_idx in spiro_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() != 'C':
            continue

        # Check if the carbon is connected to two oxygens
        oxygen_count = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() == 'O')
        if oxygen_count == 2:
            return True, "Spiroketal found"

    return False, "No spiroketal found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:72600',
                          'name': 'spiroketal',
                          'definition': 'A cyclic ketal in which the ketal '
                                        'carbon is the only common atom of two '
                                        'rings.',
                          'parents': ['CHEBI:37948', 'CHEBI:59779']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[02:32:27] SMILES Parse Error: syntax error while parsing: '
             'CC[C@@H]1CC[C@@H]2O[C@@]3(O[C@@H](C[C@@H](C)O)[C@@H](C)CC3=O)[C@@H](C)[C@H](OC(=O)\\C=C\\[C@H](C)[C@@H](O)[C@H](C)C(=O)[C@H](C)[C@@H](O)[C@H](C)C(=O)[C@](C)(O)[C@@H](O)[C@H](C)C\\C=C\\C=C\x01)[C@@H]2C\n'
             '[02:32:27] SMILES Parse Error: Failed parsing SMILES '
             "'CC[C@@H]1CC[C@@H]2O[C@@]3(O[C@@H](C[C@@H](C)O)[C@@H](C)CC3=O)[C@@H](C)[C@H](OC(=O)\\C=C\\[C@H](C)[C@@H](O)[C@H](C)C(=O)[C@H](C)[C@@H](O)[C@H](C)C(=O)[C@](C)(O)[C@@H](O)[C@H](C)C\\C=C\\C=C\x01)[C@@H]2C' "
             'for input: '
             "'CC[C@@H]1CC[C@@H]2O[C@@]3(O[C@@H](C[C@@H](C)O)[C@@H](C)CC3=O)[C@@H](C)[C@H](OC(=O)\\C=C\\[C@H](C)[C@@H](O)[C@H](C)C(=O)[C@H](C)[C@@H](O)[C@H](C)C(=O)[C@](C)(O)[C@@H](O)[C@H](C)C\\C=C\\C=C\x01)[C@@H]2C'\n",
    'stdout': '',
    'num_true_positives': 11,
    'num_false_positives': 3,
    'num_true_negatives': 9,
    'num_false_negatives': 1,
    'precision': 0.7857142857142857,
    'recall': 0.9166666666666666,
    'f1': 0.8461538461538461,
    'accuracy': None}