"""
Classifies: CHEBI:48847 heterocyclic fatty acid
"""
from rdkit import Chem

def is_heterocyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a heterocyclic fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a heterocyclic fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of carboxylic acid group (-COOH)
    carboxylic_acid = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:  # Oxygen
                    for neighbor2 in neighbor.GetNeighbors():
                        if neighbor2.GetAtomicNum() == 8 and neighbor2.GetIdx() != atom.GetIdx():
                            carboxylic_acid = True
                            break
                    if carboxylic_acid:
                        break
        if carboxylic_acid:
            break
    if not carboxylic_acid:
        return False, "No carboxylic acid group found"

    # Check for presence of ring with at least one heteroatom
    ring_info = mol.GetRingInfo()
    heterocyclic_ring = False
    for ring in ring_info.AtomRings():
        has_heteroatom = any(mol.GetAtomWithIdx(idx).GetAtomicNum() not in [6, 1] for idx in ring)
        if has_heteroatom:
            heterocyclic_ring = True
            break
    if not heterocyclic_ring:
        return False, "No heterocyclic ring found"

    return True, "Heterocyclic fatty acid"

# Example usage
smiles_list = [
    "O1C(CCCCCCCC(O)=O)=CC=C1CCCCC",
    "O=C1OC(=C)C(C1=C(O)/C(=C/CC/C=C/C=C/[C@H](CC(=O)O)C)/C)=O",
    "O=C(C=1OC(C(O)CCCC)=CC1)C(O)C",
    "C1(C(C/C=C\CCCCCO)O1)CCCCCCCC(=O)O",
    "O=C(CCC/C=C\C/C=C\CC(/C=C/[C@@H]1[C@H](CCCCC)O1)O)O",
    "CCC1OC1C\C=C/C\C=C/C\C=C/C\C=C/CCCC(O)=O",
    "OC(=O)CCCCC1CCSS1",
    "O=C(CCC/C=C\C/C=C\C/C=C\[C@@H]([C@H]1[C@H](CCCCC)O1)O)O",
    "O1C(CCCCCCCC(O)=O)=C(C(=C1CCC)C)C",
    "CCC1C(C/C=C\C/C=C\CCCCCCCC(O)=O)O1",
    "O1C(CC/C=C\CCCCCCCC/C=C\CCCCCC)=CC=C1CC(O)=O",
    "O1C(CCCCCCCCC(O)=O)=C(C(=C1CCCC)C)C",
    "O=C(C=1OC(CCCC(=O)C)=CC1)[C@@H](O)C",
    "CCCCC\C=C/C\C=C/CC1OC1C\C=C/CCCC(O)=O",
    "O=C1C=2OC(/C=C/C)=CC2C[C@H]([C@]1(OC(=O)C3=C(O)C=C(O)C=C3C)C)O",
    "C(CCC)C/C=C\C[C@@H]1[C@H](C/C=C\C/C=C\CCCC(O)=O)O1",
    "O1C(CCCCC(O)=O)=C(C(=C1/C=C/C(O)=O)C)C",
    "[H]C(\C=C/CCCCCCCC(O)=O)=C1O[C@H]1C\C=C/CC",
    "O=C(O)C(C(OC(=O)CC)C(C1OC2(OC(C3OC(C4OC(C5OC(O)(C(C)CC5C)CO)CC4C)(C)CC3)(C)CC2)CC(C1C)O)C)C",
    "O1C(=CC=C1C(OC)CCC=C)C(O)C(O)C",
    "O=C1C(=C2[C@@H](OC(C)(C)OC2)[C@@H]3[C@]1(O3)C/C=C(/C(=O)O)\C)/C=C/CCCCC",
    "OC(CCCCCC1C(/C=C/C=C/C=C\C=C\C(C/C=C\CC)O)O1)=O",
    "OC1(C2CCC(N(C2)C)C1)CC(CC(O)=O)C"
]

for smiles in smiles_list:
    result, reason = is_heterocyclic_fatty_acid(smiles)
    print(f"SMILES: {smiles} -> {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:48847',
                          'name': 'heterocyclic fatty acid',
                          'definition': 'Any fatty acid containing a ring '
                                        'composed of atoms including at least '
                                        'one heteroatom.',
                          'parents': ['CHEBI:5686', 'CHEBI:59238']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: O1C(CCCCCCCC(O)=O)=CC=C1CCCCC -> False, Reason: No '
              'carboxylic acid group found\n'
              'SMILES: '
              'O=C1OC(=C)C(C1=C(O)/C(=C/CC/C=C/C=C/[C@H](CC(=O)O)C)/C)=O -> '
              'False, Reason: No carboxylic acid group found\n'
              'SMILES: O=C(C=1OC(C(O)CCCC)=CC1)C(O)C -> False, Reason: No '
              'carboxylic acid group found\n'
              'SMILES: C1(C(C/C=C\\CCCCCO)O1)CCCCCCCC(=O)O -> False, Reason: '
              'No carboxylic acid group found\n'
              'SMILES: O=C(CCC/C=C\\C/C=C\\CC(/C=C/[C@@H]1[C@H](CCCCC)O1)O)O '
              '-> False, Reason: No carboxylic acid group found\n'
              'SMILES: CCC1OC1C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(O)=O -> False, '
              'Reason: No carboxylic acid group found\n'
              'SMILES: OC(=O)CCCCC1CCSS1 -> False, Reason: No carboxylic acid '
              'group found\n'
              'SMILES: '
              'O=C(CCC/C=C\\C/C=C\\C/C=C\\[C@@H]([C@H]1[C@H](CCCCC)O1)O)O -> '
              'False, Reason: No carboxylic acid group found\n'
              'SMILES: O1C(CCCCCCCC(O)=O)=C(C(=C1CCC)C)C -> False, Reason: No '
              'carboxylic acid group found\n'
              'SMILES: CCC1C(C/C=C\\C/C=C\\CCCCCCCC(O)=O)O1 -> False, Reason: '
              'No carboxylic acid group found\n'
              'SMILES: O1C(CC/C=C\\CCCCCCCC/C=C\\CCCCCC)=CC=C1CC(O)=O -> '
              'False, Reason: No carboxylic acid group found\n'
              'SMILES: O1C(CCCCCCCCC(O)=O)=C(C(=C1CCCC)C)C -> False, Reason: '
              'No carboxylic acid group found\n'
              'SMILES: O=C(C=1OC(CCCC(=O)C)=CC1)[C@@H](O)C -> False, Reason: '
              'No carboxylic acid group found\n'
              'SMILES: CCCCC\\C=C/C\\C=C/CC1OC1C\\C=C/CCCC(O)=O -> False, '
              'Reason: No carboxylic acid group found\n'
              'SMILES: '
              'O=C1C=2OC(/C=C/C)=CC2C[C@H]([C@]1(OC(=O)C3=C(O)C=C(O)C=C3C)C)O '
              '-> False, Reason: No carboxylic acid group found\n'
              'SMILES: C(CCC)C/C=C\\C[C@@H]1[C@H](C/C=C\\C/C=C\\CCCC(O)=O)O1 '
              '-> False, Reason: No carboxylic acid group found\n'
              'SMILES: O1C(CCCCC(O)=O)=C(C(=C1/C=C/C(O)=O)C)C -> False, '
              'Reason: No carboxylic acid group found\n'
              'SMILES: [H]C(\\C=C/CCCCCCCC(O)=O)=C1O[C@H]1C\\C=C/CC -> False, '
              'Reason: No carboxylic acid group found\n'
              'SMILES: '
              'O=C(O)C(C(OC(=O)CC)C(C1OC2(OC(C3OC(C4OC(C5OC(O)(C(C)CC5C)CO)CC4C)(C)CC3)(C)CC2)CC(C1C)O)C)C '
              '-> False, Reason: No carboxylic acid group found\n'
              'SMILES: O1C(=CC=C1C(OC)CCC=C)C(O)C(O)C -> False, Reason: No '
              'carboxylic acid group found\n'
              'SMILES: '
              'O=C1C(=C2[C@@H](OC(C)(C)OC2)[C@@H]3[C@]1(O3)C/C=C(/C(=O)O)\\C)/C=C/CCCCC '
              '-> False, Reason: No carboxylic acid group found\n'
              'SMILES: OC(CCCCCC1C(/C=C/C=C/C=C\\C=C\\C(C/C=C\\CC)O)O1)=O -> '
              'False, Reason: No carboxylic acid group found\n'
              'SMILES: OC1(C2CCC(N(C2)C)C1)CC(CC(O)=O)C -> False, Reason: No '
              'carboxylic acid group found\n',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 23,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}