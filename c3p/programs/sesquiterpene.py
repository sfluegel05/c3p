"""
Classifies: CHEBI:35189 sesquiterpene
"""
from rdkit import Chem

def is_sesquiterpene(smiles: str):
    """
    Determines if a molecule is a sesquiterpene (C15 terpene).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')

    # Sesquiterpenes should have exactly 15 carbon atoms
    if num_carbons != 15:
        return False, f"Number of carbon atoms is {num_carbons}, not 15"

    return True, "Molecule is a sesquiterpene"

# Examples for testing
examples = [
    "C[C@@H]1CC[C@@]2(C)C(C)=C[C@]3(C)CCC[C@]123",
    "C\\C1=C/CC\\C(C)=C\\CC(CC1)=C(C)C",
    "[H][C@@]12CCC(=C)[C@]1([H])[C@@]1([H])C(C)(C)CC[C@@]1(C)C2",
    "C[C@@H]1CC[C@@H](C(C)=C)[C@]11CCC(C)=CC1",
    "C\\C=C(/C)CC\\C=C(/C)CCC=C(C)C",
    "CC(C)=CC\\C=C(/C)C1CCC(C)=CC1",
    "C1C[C@H](C=C2[C@@]([C@@H]1C)(CCC2C)[H])C(=C)C",
    "CC(C)=CCC\\C(C)=C/CCC(=C)C=C",
    "[H]C(CCC(=C)C=C)=C(C)CCC=C(C)C",
    "CC1(C)CCCC2C(C)(C)[C@H]3CC[C@@]12C3",
    "CC(C)=CCCC1(C)C2CC=C(C)C1C2",
    "[H][C@]12CCC(C)=C1[C@@]1([H])[C@@]([H])(CC[C@H]2C)C1(C)C",
    "[H][C@@]1(CCC(C)=CC1)C(=C)CCC=C(C)C",
    "C1C(C2[C@@]3([C@@H](CC[C@@]2([C@@]3(C1)[H])C)C(C)C)[H])=C",
    "[H][C@@]12CCCCC1C1CC[C@H]3CC[C@@]1(C3)C2",
    "CC(C)=CCC[C@@]1(C)[C@H]2CC=C(C)[C@@H]1C2",
    "[C@@]12(C([C@]1(CC[C@@H]([C@]3(C2=C(CC3)C)[H])C)[H])(C)C)[H]",
    "[H][C@@]12CCC(C)=C[C@@]1([H])C(C)(C)CCC=C2C",
    "[H][C@@]1(CC=C(C)C=C1)[C@@H](C)CCC=C(C)C",
    "[H][C@]12CC[C@]3([H])[C@@]1([H])C(C)(C)CCC[C@]3(C)C2=C",
    "[H][C@]1(CC[C@@]2(C)CCC=C(C)[C@]2([H])C1)C(C)=C",
    "CC(C)=CCC\\C(C)=C1/CCC(C)=CC1"
]

for example in examples:
    result, reason = is_sesquiterpene(example)
    print(f"SMILES: {example} -> {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35189',
                          'name': 'sesquiterpene',
                          'definition': 'A C15 terpene.',
                          'parents': ['CHEBI:35186']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[23:23:25] SMILES Parse Error: unclosed ring for input: '
             "'CC(C)C1=C/C=C(C)\\CC\\C=C(C)\\CC\\C=C(C)\\CC\x01'\n",
    'stdout': 'SMILES: C[C@@H]1CC[C@@]2(C)C(C)=C[C@]3(C)CCC[C@]123 -> True, '
              'Reason: Molecule is a sesquiterpene\n'
              'SMILES: C\\C1=C/CC\\C(C)=C\\CC(CC1)=C(C)C -> True, Reason: '
              'Molecule is a sesquiterpene\n'
              'SMILES: '
              '[H][C@@]12CCC(=C)[C@]1([H])[C@@]1([H])C(C)(C)CC[C@@]1(C)C2 -> '
              'True, Reason: Molecule is a sesquiterpene\n'
              'SMILES: C[C@@H]1CC[C@@H](C(C)=C)[C@]11CCC(C)=CC1 -> True, '
              'Reason: Molecule is a sesquiterpene\n'
              'SMILES: C\\C=C(/C)CC\\C=C(/C)CCC=C(C)C -> True, Reason: '
              'Molecule is a sesquiterpene\n'
              'SMILES: CC(C)=CC\\C=C(/C)C1CCC(C)=CC1 -> True, Reason: Molecule '
              'is a sesquiterpene\n'
              'SMILES: C1C[C@H](C=C2[C@@]([C@@H]1C)(CCC2C)[H])C(=C)C -> True, '
              'Reason: Molecule is a sesquiterpene\n'
              'SMILES: CC(C)=CCC\\C(C)=C/CCC(=C)C=C -> True, Reason: Molecule '
              'is a sesquiterpene\n'
              'SMILES: [H]C(CCC(=C)C=C)=C(C)CCC=C(C)C -> True, Reason: '
              'Molecule is a sesquiterpene\n'
              'SMILES: CC1(C)CCCC2C(C)(C)[C@H]3CC[C@@]12C3 -> True, Reason: '
              'Molecule is a sesquiterpene\n'
              'SMILES: CC(C)=CCCC1(C)C2CC=C(C)C1C2 -> True, Reason: Molecule '
              'is a sesquiterpene\n'
              'SMILES: '
              '[H][C@]12CCC(C)=C1[C@@]1([H])[C@@]([H])(CC[C@H]2C)C1(C)C -> '
              'True, Reason: Molecule is a sesquiterpene\n'
              'SMILES: [H][C@@]1(CCC(C)=CC1)C(=C)CCC=C(C)C -> True, Reason: '
              'Molecule is a sesquiterpene\n'
              'SMILES: '
              'C1C(C2[C@@]3([C@@H](CC[C@@]2([C@@]3(C1)[H])C)C(C)C)[H])=C -> '
              'True, Reason: Molecule is a sesquiterpene\n'
              'SMILES: [H][C@@]12CCCCC1C1CC[C@H]3CC[C@@]1(C3)C2 -> True, '
              'Reason: Molecule is a sesquiterpene\n'
              'SMILES: CC(C)=CCC[C@@]1(C)[C@H]2CC=C(C)[C@@H]1C2 -> True, '
              'Reason: Molecule is a sesquiterpene\n'
              'SMILES: '
              '[C@@]12(C([C@]1(CC[C@@H]([C@]3(C2=C(CC3)C)[H])C)[H])(C)C)[H] -> '
              'True, Reason: Molecule is a sesquiterpene\n'
              'SMILES: [H][C@@]12CCC(C)=C[C@@]1([H])C(C)(C)CCC=C2C -> True, '
              'Reason: Molecule is a sesquiterpene\n'
              'SMILES: [H][C@@]1(CC=C(C)C=C1)[C@@H](C)CCC=C(C)C -> True, '
              'Reason: Molecule is a sesquiterpene\n'
              'SMILES: [H][C@]12CC[C@]3([H])[C@@]1([H])C(C)(C)CCC[C@]3(C)C2=C '
              '-> True, Reason: Molecule is a sesquiterpene\n'
              'SMILES: [H][C@]1(CC[C@@]2(C)CCC=C(C)[C@]2([H])C1)C(C)=C -> '
              'True, Reason: Molecule is a sesquiterpene\n'
              'SMILES: CC(C)=CCC\\C(C)=C1/CCC(C)=CC1 -> True, Reason: Molecule '
              'is a sesquiterpene\n',
    'num_true_positives': 22,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}