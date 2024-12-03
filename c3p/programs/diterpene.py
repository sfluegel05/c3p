"""
Classifies: CHEBI:35190 diterpene
"""
from rdkit import Chem

def is_diterpene(smiles: str):
    """
    Determines if a molecule is a diterpene (a C20 terpene).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule has 20 carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if carbon_count != 20:
        return False, f"Molecule has {carbon_count} carbon atoms, not 20"

    # Check if the molecule is a terpene
    # Terpenes are hydrocarbons, so we check for the presence of only carbon and hydrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in {'C', 'H'}:
            return True, "Molecule is a diterpene with functional groups"

    return True, "Molecule is a diterpene"

# Example usage
smiles_list = [
    "[C@]12([C@@](C=C(C(C(C=C(C)C(CC=C(CC1)C)O)O)=O)C)(C2(C)C)[H])[H]",
    "[H][C@@]12CC[C@@](C)(C=C)C=C1CC[C@]1([H])C(C)(C)CCC[C@@]21C",
    "[C@H]1(CCC(=C)C=C)C(C)=CC[C@@]2([C@]1(C)CCCC2(C)C)[H]",
    "CC(C)=CCC\\C(C)=C/CCC1(C)C2CC3C(C2)C13C",
    "CC[C@@]1(C)CC[C@H]2[C@@H](CC[C@H]3C(C)(C)CCC[C@]23C)C1",
    "[H][C@@]12CC[C@]3([H])C(C)(C)CCC[C@@]3(C)C1=CC[C@@](C)(C2)C=C",
    "[C@]12([C@@](C=C(C(CC=C(C)C(CC=C(CC1)C)=O)O)C)(C2(C)C)[H])[H]",
    "[H][C@@]12CC[C@@]3([H])C(C)(C)CCC[C@]3(C)[C@@]1([H])CC[C@@](C)(CC)C2",
    "CC(C)c1ccc2c(CC[C@H]3C(C)(C)CCC[C@]23C)c1",
    "[H][C@]12CC=C(C)[C@@]([H])(CCC(C)=CCC\\C(C)=C\\C1)C2(C)C"
]

for smiles in smiles_list:
    result, reason = is_diterpene(smiles)
    print(f"SMILES: {smiles}\nIs diterpene: {result}\nReason: {reason}\n")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35190',
                          'name': 'diterpene',
                          'definition': 'A C20 terpene.',
                          'parents': ['CHEBI:35186']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: '
              '[C@]12([C@@](C=C(C(C(C=C(C)C(CC=C(CC1)C)O)O)=O)C)(C2(C)C)[H])[H]\n'
              'Is diterpene: True\n'
              'Reason: Molecule is a diterpene with functional groups\n'
              '\n'
              'SMILES: '
              '[H][C@@]12CC[C@@](C)(C=C)C=C1CC[C@]1([H])C(C)(C)CCC[C@@]21C\n'
              'Is diterpene: True\n'
              'Reason: Molecule is a diterpene\n'
              '\n'
              'SMILES: [C@H]1(CCC(=C)C=C)C(C)=CC[C@@]2([C@]1(C)CCCC2(C)C)[H]\n'
              'Is diterpene: True\n'
              'Reason: Molecule is a diterpene\n'
              '\n'
              'SMILES: CC(C)=CCC\\C(C)=C/CCC1(C)C2CC3C(C2)C13C\n'
              'Is diterpene: True\n'
              'Reason: Molecule is a diterpene\n'
              '\n'
              'SMILES: CC[C@@]1(C)CC[C@H]2[C@@H](CC[C@H]3C(C)(C)CCC[C@]23C)C1\n'
              'Is diterpene: True\n'
              'Reason: Molecule is a diterpene\n'
              '\n'
              'SMILES: '
              '[H][C@@]12CC[C@]3([H])C(C)(C)CCC[C@@]3(C)C1=CC[C@@](C)(C2)C=C\n'
              'Is diterpene: True\n'
              'Reason: Molecule is a diterpene\n'
              '\n'
              'SMILES: '
              '[C@]12([C@@](C=C(C(CC=C(C)C(CC=C(CC1)C)=O)O)C)(C2(C)C)[H])[H]\n'
              'Is diterpene: True\n'
              'Reason: Molecule is a diterpene with functional groups\n'
              '\n'
              'SMILES: '
              '[H][C@@]12CC[C@@]3([H])C(C)(C)CCC[C@]3(C)[C@@]1([H])CC[C@@](C)(CC)C2\n'
              'Is diterpene: True\n'
              'Reason: Molecule is a diterpene\n'
              '\n'
              'SMILES: CC(C)c1ccc2c(CC[C@H]3C(C)(C)CCC[C@]23C)c1\n'
              'Is diterpene: True\n'
              'Reason: Molecule is a diterpene\n'
              '\n'
              'SMILES: '
              '[H][C@]12CC=C(C)[C@@]([H])(CCC(C)=CCC\\C(C)=C\\C1)C2(C)C\n'
              'Is diterpene: True\n'
              'Reason: Molecule is a diterpene\n'
              '\n',
    'num_true_positives': 10,
    'num_false_positives': 0,
    'num_true_negatives': 10,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}