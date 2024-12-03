"""
Classifies: CHEBI:26872 terpene ketone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_terpene_ketone(smiles: str):
    """
    Determines if a molecule is a terpene ketone (a terpenoid which contains a keto group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a terpene ketone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule contains a keto group (C=O)
    keto_group = Chem.MolFromSmarts('C=O')
    if not mol.HasSubstructMatch(keto_group):
        return False, "No keto group found"

    # Check if molecule is a terpenoid (contains isoprene units: C5H8)
    isoprene_unit = Chem.MolFromSmarts('C(C)(C)=C')
    if not mol.HasSubstructMatch(isoprene_unit):
        return False, "No isoprene units found, not a terpenoid"

    return True, "Molecule is a terpene ketone"

# Example usage
smiles_examples = [
    "[C@@]12(CC[C@]3([C@](C[C@H](C([C@@H]3C)=O)O)([H])[C@@]1(CC[C@@]4([C@@]2(CC[C@@]5([C@]4(C[C@@](CC5)(C(O)=O)C)[H])C)C)C)C)C)[H]",
    "CC(=O)O[C@@H]1C[C@H]2C(C)(C)C(=O)C=C[C@]2(C)[C@H]2CC[C@]3(C)C(=CC(=O)O[C@]3(O)c3ccoc3)[C@]12C",
    "CC(C)CCC[C@H](C)[C@@H]1CC[C@]2(C)[C@@H]3CC=C4C(C)(C)C(=O)C[C@@H](O)[C@]4(C)[C@H]3CC[C@@]12C",
    "CC(C)[C@@H]1C[C@@H](O)[C@H]2[C@@]1(CO)CC[C@@]1(C)[C@@H]3[C@@H](O)C[C@H]4C(C)(C)C(=O)CC[C@]4(C)C3=CC[C@]21C",
    "CC1(C)[C@@H]2CC[C@@]1(C)C(=O)C2",
    "CC(=O)CC\C=C(/C)CC\C=C(/C)CCC=C(C)C"
]

for smiles in smiles_examples:
    result, reason = is_terpene_ketone(smiles)
    print(f"SMILES: {smiles}\nResult: {result}\nReason: {reason}\n")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26872',
                          'name': 'terpene ketone',
                          'definition': 'Any terpenoid which contains a keto '
                                        'group.',
                          'parents': ['CHEBI:17087', 'CHEBI:26873']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: '
              '[C@@]12(CC[C@]3([C@](C[C@H](C([C@@H]3C)=O)O)([H])[C@@]1(CC[C@@]4([C@@]2(CC[C@@]5([C@]4(C[C@@](CC5)(C(O)=O)C)[H])C)C)C)C)C)[H]\n'
              'Result: False\n'
              'Reason: No isoprene units found, not a terpenoid\n'
              '\n'
              'SMILES: '
              'CC(=O)O[C@@H]1C[C@H]2C(C)(C)C(=O)C=C[C@]2(C)[C@H]2CC[C@]3(C)C(=CC(=O)O[C@]3(O)c3ccoc3)[C@]12C\n'
              'Result: True\n'
              'Reason: Molecule is a terpene ketone\n'
              '\n'
              'SMILES: '
              'CC(C)CCC[C@H](C)[C@@H]1CC[C@]2(C)[C@@H]3CC=C4C(C)(C)C(=O)C[C@@H](O)[C@]4(C)[C@H]3CC[C@@]12C\n'
              'Result: True\n'
              'Reason: Molecule is a terpene ketone\n'
              '\n'
              'SMILES: '
              'CC(C)[C@@H]1C[C@@H](O)[C@H]2[C@@]1(CO)CC[C@@]1(C)[C@@H]3[C@@H](O)C[C@H]4C(C)(C)C(=O)CC[C@]4(C)C3=CC[C@]21C\n'
              'Result: True\n'
              'Reason: Molecule is a terpene ketone\n'
              '\n'
              'SMILES: CC1(C)[C@@H]2CC[C@@]1(C)C(=O)C2\n'
              'Result: False\n'
              'Reason: No isoprene units found, not a terpenoid\n'
              '\n'
              'SMILES: CC(=O)CC\\C=C(/C)CC\\C=C(/C)CCC=C(C)C\n'
              'Result: True\n'
              'Reason: Molecule is a terpene ketone\n'
              '\n',
    'num_true_positives': 15,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 5,
    'precision': 0.9375,
    'recall': 0.75,
    'f1': 0.8333333333333334,
    'accuracy': None}