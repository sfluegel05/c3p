"""
Classifies: CHEBI:33654 alicyclic compound
"""
from rdkit import Chem

def is_alicyclic_compound(smiles: str):
    """
    Determines if a molecule is an alicyclic compound (aliphatic compound with a carbocyclic ring structure,
    which may be saturated or unsaturated, but not a benzenoid or other aromatic system).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alicyclic compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check if there is at least one ring
    if not rings.AtomRings():
        return False, "No rings found"

    # Check for aromaticity
    for ring in rings.AtomRings():
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if any(atom.GetIsAromatic() for atom in atoms):
            return False, "Contains aromatic ring"

    # Check if all atoms in the ring are carbons
    for ring in rings.AtomRings():
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if not all(atom.GetSymbol() == 'C' for atom in atoms):
            return False, "Ring contains non-carbon atoms"

    return True, "Alicyclic compound"

# Example usage
smiles_examples = [
    "OC1C=CC=C(C1O)C(O)=O",  # 2,3-dihydroxy-2,3-dihydrobenzoic acid
    "C=CCC1CCCCC1",  # Prop-2-enylcyclohexane
    "[H][C@@]12CC[C@H](C(C)C)[C@@]1(C)C\\C=C(C)/C1=CC[C@@H](C)[C@@]1([H])C2",  # collinodiene
    "CC1=CC=CC(O)C1(O)C(O)=O",  # 1,6-dihydroxy-2-methylcyclohexa-2,4-dienecarboxylic acid
    "[H][C@@]12\\C=C\\C(=C)CCC\\C(C)=C\\CC\\C(C)=C\\C[C@@]1(C)CC[C@@H]2C(C)C",  # fusaproliferene
    "[H][C@]12C[C@@]3(C)CC[C@]4(C)CC[C@H](C(C)C)[C@@]4([H])[C@]3([H])[C@@]1([H])C[C@]1(C)CC[C@H]2C1=C",  # subrutilene F
    "OC=1CCCC(C1)(C)C",  # 3,3-dimethyl-1-cyclohexen-1-ol
    "[H][C@]12C[C@@H](O)C(CO)=C[C@H](O)[C@]1(C)CC[C@@]1(C)CCC(C(C)C)=C21",  # cyathatriol
    "C1C=CC=C1",  # cyclopentadiene
    "CC(C)C1=CC=C(C)CC1",  # alpha-terpinene
    "CC=1CCCC(C1)=O",  # 3-Methyl-2-cyclohexen-1-one
    "O[C@H]1C=CC=C[C@]1(O)C([O-])=O",  # (1R,6S)-1,6-dihydroxycyclohexa-2,4-diene-1-carboxylate
    "C1(CCCCC1)(C)C",  # 1,1-dimethyl-Cyclohexane
    "[H][C@@]12C=C(C)CC[C@@]11[C@H](C)CCC1=C(C)CC[C@@H]2C(C)C",  # bonnadiene
    "CC(C)=CCCC(=C)C1CCC=C(CCC=C(C)C)C1",  # gamma-camphorene
    "C1C[C@@]2([C@@](C=3[C@@]1(CCC3C(C)C)C)(C[C@H](C(=C[C@@H]2O[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)CO)OC(C)=O)[H])C",  # erinacine Q2
    "C=C1C=C(C)CC(C)(C)C1",  # 3-Methylene-1,5,5-trimethylcyclohexene
    "C[C@H]1CC\\C=C(C)/CC2=C1CC(C)(C)C2",  # (3S)-(+)-asterisca-2(9),6-diene
    "C1(CCCC=C1)CC",  # 3-ethyl-Cyclohexene
    "[H][C@@]1(CC=C(C)C=C1)[C@H](C)CCC=C(C)C",  # 7-epi-zingiberene
    "[C@@]12(C(CCC(=C[C@@]1([C@@H](CC2)C)[H])C(C)C)=C)[H]",  # (1R,4R,5S)-(-)-guaia-6,10(14)-diene
    "C1(CCCCC1)CCCCCCCCCCC",  # undecyl-Cyclohexane
    "C1(CC(CC1)CCCCCCCCC)CCC(CCCCCCCCCCCC)C",  # 1-(3-methylpentadecyl)-3-nonylcyclopentane
    "[C@@H]1([C@@H](CCCC1)C)C",  # trans-1,2-dimethylcyclohexane
    "[H][C@@]1(CCC(=C)C=C1)C(C)C",  # (-)-beta-phellandrene
    "[H][C@]12C[C@]3(C)CCC[C@]3(C)C1=C[C@@]1(C)CC[C@]3(C)CC[C@H](C(C)C)[C@@]3([H])[C@]21[H]",  # subrutilene D
    "CC(C1CCC([C@]12CCC3(C4CCC(C4C23)C)C)=C)C"  # spiroluchuene A
]

for smiles in smiles_examples:
    result, reason = is_alicyclic_compound(smiles)
    print(f"SMILES: {smiles}, Result: {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33654',
                          'name': 'alicyclic compound',
                          'definition': 'An aliphatic compound having a '
                                        'carbocyclic ring structure which may '
                                        'be saturated or unsaturated, but may '
                                        'not be a benzenoid or other aromatic '
                                        'system.',
                          'parents': ['CHEBI:33598', 'CHEBI:33653']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: OC1C=CC=C(C1O)C(O)=O, Result: True, Reason: Alicyclic '
              'compound\n'
              'SMILES: C=CCC1CCCCC1, Result: True, Reason: Alicyclic compound\n'
              'SMILES: '
              '[H][C@@]12CC[C@H](C(C)C)[C@@]1(C)C\\C=C(C)/C1=CC[C@@H](C)[C@@]1([H])C2, '
              'Result: True, Reason: Alicyclic compound\n'
              'SMILES: CC1=CC=CC(O)C1(O)C(O)=O, Result: True, Reason: '
              'Alicyclic compound\n'
              'SMILES: '
              '[H][C@@]12\\C=C\\C(=C)CCC\\C(C)=C\\CC\\C(C)=C\\C[C@@]1(C)CC[C@@H]2C(C)C, '
              'Result: True, Reason: Alicyclic compound\n'
              'SMILES: '
              '[H][C@]12C[C@@]3(C)CC[C@]4(C)CC[C@H](C(C)C)[C@@]4([H])[C@]3([H])[C@@]1([H])C[C@]1(C)CC[C@H]2C1=C, '
              'Result: True, Reason: Alicyclic compound\n'
              'SMILES: OC=1CCCC(C1)(C)C, Result: True, Reason: Alicyclic '
              'compound\n'
              'SMILES: '
              '[H][C@]12C[C@@H](O)C(CO)=C[C@H](O)[C@]1(C)CC[C@@]1(C)CCC(C(C)C)=C21, '
              'Result: True, Reason: Alicyclic compound\n'
              'SMILES: C1C=CC=C1, Result: True, Reason: Alicyclic compound\n'
              'SMILES: CC(C)C1=CC=C(C)CC1, Result: True, Reason: Alicyclic '
              'compound\n'
              'SMILES: CC=1CCCC(C1)=O, Result: True, Reason: Alicyclic '
              'compound\n'
              'SMILES: O[C@H]1C=CC=C[C@]1(O)C([O-])=O, Result: True, Reason: '
              'Alicyclic compound\n'
              'SMILES: C1(CCCCC1)(C)C, Result: True, Reason: Alicyclic '
              'compound\n'
              'SMILES: '
              '[H][C@@]12C=C(C)CC[C@@]11[C@H](C)CCC1=C(C)CC[C@@H]2C(C)C, '
              'Result: True, Reason: Alicyclic compound\n'
              'SMILES: CC(C)=CCCC(=C)C1CCC=C(CCC=C(C)C)C1, Result: True, '
              'Reason: Alicyclic compound\n'
              'SMILES: '
              'C1C[C@@]2([C@@](C=3[C@@]1(CCC3C(C)C)C)(C[C@H](C(=C[C@@H]2O[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)O)CO)OC(C)=O)[H])C, '
              'Result: False, Reason: Ring contains non-carbon atoms\n'
              'SMILES: C=C1C=C(C)CC(C)(C)C1, Result: True, Reason: Alicyclic '
              'compound\n'
              'SMILES: C[C@H]1CC\\C=C(C)/CC2=C1CC(C)(C)C2, Result: True, '
              'Reason: Alicyclic compound\n'
              'SMILES: C1(CCCC=C1)CC, Result: True, Reason: Alicyclic '
              'compound\n'
              'SMILES: [H][C@@]1(CC=C(C)C=C1)[C@H](C)CCC=C(C)C, Result: True, '
              'Reason: Alicyclic compound\n'
              'SMILES: [C@@]12(C(CCC(=C[C@@]1([C@@H](CC2)C)[H])C(C)C)=C)[H], '
              'Result: True, Reason: Alicyclic compound\n'
              'SMILES: C1(CCCCC1)CCCCCCCCCCC, Result: True, Reason: Alicyclic '
              'compound\n'
              'SMILES: C1(CC(CC1)CCCCCCCCC)CCC(CCCCCCCCCCCC)C, Result: True, '
              'Reason: Alicyclic compound\n'
              'SMILES: [C@@H]1([C@@H](CCCC1)C)C, Result: True, Reason: '
              'Alicyclic compound\n'
              'SMILES: [H][C@@]1(CCC(=C)C=C1)C(C)C, Result: True, Reason: '
              'Alicyclic compound\n'
              'SMILES: '
              '[H][C@]12C[C@]3(C)CCC[C@]3(C)C1=C[C@@]1(C)CC[C@]3(C)CC[C@H](C(C)C)[C@@]3([H])[C@]21[H], '
              'Result: True, Reason: Alicyclic compound\n'
              'SMILES: CC(C1CCC([C@]12CCC3(C4CCC(C4C23)C)C)=C)C, Result: True, '
              'Reason: Alicyclic compound\n',
    'num_true_positives': 26,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9629629629629629,
    'f1': 0.9811320754716981,
    'accuracy': None}