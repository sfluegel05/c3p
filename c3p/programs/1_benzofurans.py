"""
Classifies: CHEBI:38830 1-benzofurans
"""
from rdkit import Chem

def is_1_benzofurans(smiles: str):
    """
    Determines if a molecule is a 1-benzofuran or its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-benzofuran or its substituted derivatives, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMILES pattern for the 1-benzofuran core
    benzofuran_pattern = Chem.MolFromSmarts("c1ccc2occc2c1")

    if benzofuran_pattern is None:
        return False, "Invalid benzofuran pattern"

    # Check if the molecule contains the 1-benzofuran core
    if mol.HasSubstructMatch(benzofuran_pattern):
        return True, "Molecule is a 1-benzofuran or its substituted derivative"
    else:
        return False, "Molecule does not contain the 1-benzofuran core"

# Example usage
smiles_examples = [
    "O=C1O[C@](C(=O)OCC)(CC2=CC3=C(O[C@H](C3)C(O)(C)C)C=C2)C(=C1O)C4=CC=C(O)C=C4",
    "Oc1ccc(\\C=C/c2ccc3O[C@H]([C@@H](c3c2)c2cc(O)cc(O)c2)c2ccc(O)cc2)cc1",
    "[H]C(c1cc2cc(OCCOc3cc(C)ccc3N(CC(=O)OCOC(C)=O)CC(=O)OCOC(C)=O)c(cc2o1)N(CC(=O)OCOC(C)=O)CC(=O)OCOC(C)=O)=C1NC(=S)NC1=O",
    "[H][C@@]12[C@H](Oc3cc(O)cc(c13)[C@@]1([H])[C@@H](Oc3cc(O)cc(c13)[C@@]1([H])[C@@H](Oc3cc(O)cc2c13)c1ccc(O)cc1)c1ccc(O)cc1)c1ccc(O)cc1",
    "[H][C@]12CC(C)=C[C@]3([H])c4c(O)cc(cc4O[C@@](Oc4cc(O)ccc14)(c1ccc(O)c4C=CC(C)(C)Oc14)[C@@]23[H])-c1cc2ccc(O)cc2o1",
    "COC(=O)C1=CC(=O)C=C(OC)C11Oc2cc(C)cc(O)c2C1=O",
    "C=1C=C2C(=CC1O)CC(=O)O2",
    "O=C1NC(=O)C(=C1C2=CC=C(O)C=C2)CC3=CC4=C(O[C@H](C4)C(O)(C)C)C=C3",
    "O1C(O)(CC2=CC(O)=C(OC)C(O)=C2)C(=O)C=3C1=CC(O)=CC3O",
    "C=12OC(C(C1C=C(C=C2OC)/C=C/C(O)=O)COS(O)(=O)=O)C3=CC(=C(C=C3)O)OC",
    "O1C2=C(C=C(COCC3=CC=C(OCC=C=C)C=C3)C=C2)C[C@H]1C(=C)C",
    "Oc1ccc(\\C=C\\c2cc(O)cc3O[C@@H]([C@H](c23)c2cc(O)cc(O)c2)c2ccc(O)cc2)cc1",
    "OCCCc1ccc2oc(cc2c1)-c1ccc2OCOc2c1",
    "C=1C2=C(C(=C(C1)[C@H](CCCC)O)CO)O[C@@H](C2)C(C)(C)O",
    "COc1c2O\\C(=C/c3ccc(O)c(O)c3)C(=O)c2ccc1O",
    "C1=C(C=C2C(=C1OC)O[C@H]([C@@H]2CO)C3=CC=C(C(=C3)OC)O)/C=C/CO",
    "O=C(C1=C(O)C2=C(OCC2)C=C1OC)/C=C/C",
    "Oc1ccc(cc1)[C@H]1O[C@H]([C@@H]([C@@H]1c1cc(O)cc(O)c1)c1cc2[C@H]([C@@H](Oc2cc1O)c1ccc(O)cc1)c1cc2[C@H]([C@@H](Oc2cc1O)c1ccc(O)cc1)c1cc(O)cc(O)c1)c1ccc(O)cc1",
    "Oc1ccc(cc1)[C@@H]1Oc2ccc(\\C=C\\c3cc(O)cc(O)c3)cc2[C@H]1c1cc(O)cc(O)c1",
    "O1[C@@H]([C@H](C2=C1C(OC)=CC(=C2)C(=O)[H])C)C3=CC(OC)=C(OC)C=C3",
    "COc1cc(O)cc2OC(O)(Cc3ccc(O)cc3)C(=O)c12",
    "CC1(CC2=C(O1)C(=CC=C2)OCC(CNCC3=CC=CC=N3)O)C",
    "C\\C=C\\c1c(O)c(O)cc2OC(=O)[C@](C)(NC(=O)\\C=C\\C(O)=O)c12",
    "COc1ccc2cc(oc2c1)-c1ccc(O)cc1OC"
]

for smiles in smiles_examples:
    result, reason = is_1_benzofurans(smiles)
    print(f"SMILES: {smiles} -> {result}, {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38830',
                          'name': '1-benzofurans',
                          'definition': 'A member of the class of benzofurans '
                                        'consisting of a 1-benzofuran skeleton '
                                        'and its substituted derivatives '
                                        'thereof.',
                          'parents': ['CHEBI:35259']},
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
              'O=C1O[C@](C(=O)OCC)(CC2=CC3=C(O[C@H](C3)C(O)(C)C)C=C2)C(=C1O)C4=CC=C(O)C=C4 '
              '-> False, Molecule does not contain the 1-benzofuran core\n'
              'SMILES: '
              'Oc1ccc(\\C=C/c2ccc3O[C@H]([C@@H](c3c2)c2cc(O)cc(O)c2)c2ccc(O)cc2)cc1 '
              '-> False, Molecule does not contain the 1-benzofuran core\n'
              'SMILES: '
              '[H]C(c1cc2cc(OCCOc3cc(C)ccc3N(CC(=O)OCOC(C)=O)CC(=O)OCOC(C)=O)c(cc2o1)N(CC(=O)OCOC(C)=O)CC(=O)OCOC(C)=O)=C1NC(=S)NC1=O '
              '-> True, Molecule is a 1-benzofuran or its substituted '
              'derivative\n'
              'SMILES: '
              '[H][C@@]12[C@H](Oc3cc(O)cc(c13)[C@@]1([H])[C@@H](Oc3cc(O)cc(c13)[C@@]1([H])[C@@H](Oc3cc(O)cc2c13)c1ccc(O)cc1)c1ccc(O)cc1)c1ccc(O)cc1 '
              '-> False, Molecule does not contain the 1-benzofuran core\n'
              'SMILES: '
              '[H][C@]12CC(C)=C[C@]3([H])c4c(O)cc(cc4O[C@@](Oc4cc(O)ccc14)(c1ccc(O)c4C=CC(C)(C)Oc14)[C@@]23[H])-c1cc2ccc(O)cc2o1 '
              '-> True, Molecule is a 1-benzofuran or its substituted '
              'derivative\n'
              'SMILES: COC(=O)C1=CC(=O)C=C(OC)C11Oc2cc(C)cc(O)c2C1=O -> False, '
              'Molecule does not contain the 1-benzofuran core\n'
              'SMILES: C=1C=C2C(=CC1O)CC(=O)O2 -> False, Molecule does not '
              'contain the 1-benzofuran core\n'
              'SMILES: '
              'O=C1NC(=O)C(=C1C2=CC=C(O)C=C2)CC3=CC4=C(O[C@H](C4)C(O)(C)C)C=C3 '
              '-> False, Molecule does not contain the 1-benzofuran core\n'
              'SMILES: O1C(O)(CC2=CC(O)=C(OC)C(O)=C2)C(=O)C=3C1=CC(O)=CC3O -> '
              'False, Molecule does not contain the 1-benzofuran core\n'
              'SMILES: '
              'C=12OC(C(C1C=C(C=C2OC)/C=C/C(O)=O)COS(O)(=O)=O)C3=CC(=C(C=C3)O)OC '
              '-> False, Molecule does not contain the 1-benzofuran core\n'
              'SMILES: O1C2=C(C=C(COCC3=CC=C(OCC=C=C)C=C3)C=C2)C[C@H]1C(=C)C '
              '-> False, Molecule does not contain the 1-benzofuran core\n'
              'SMILES: '
              'Oc1ccc(\\C=C\\c2cc(O)cc3O[C@@H]([C@H](c23)c2cc(O)cc(O)c2)c2ccc(O)cc2)cc1 '
              '-> False, Molecule does not contain the 1-benzofuran core\n'
              'SMILES: OCCCc1ccc2oc(cc2c1)-c1ccc2OCOc2c1 -> True, Molecule is '
              'a 1-benzofuran or its substituted derivative\n'
              'SMILES: C=1C2=C(C(=C(C1)[C@H](CCCC)O)CO)O[C@@H](C2)C(C)(C)O -> '
              'False, Molecule does not contain the 1-benzofuran core\n'
              'SMILES: COc1c2O\\C(=C/c3ccc(O)c(O)c3)C(=O)c2ccc1O -> False, '
              'Molecule does not contain the 1-benzofuran core\n'
              'SMILES: '
              'C1=C(C=C2C(=C1OC)O[C@H]([C@@H]2CO)C3=CC=C(C(=C3)OC)O)/C=C/CO -> '
              'False, Molecule does not contain the 1-benzofuran core\n'
              'SMILES: O=C(C1=C(O)C2=C(OCC2)C=C1OC)/C=C/C -> False, Molecule '
              'does not contain the 1-benzofuran core\n'
              'SMILES: '
              'Oc1ccc(cc1)[C@H]1O[C@H]([C@@H]([C@@H]1c1cc(O)cc(O)c1)c1cc2[C@H]([C@@H](Oc2cc1O)c1ccc(O)cc1)c1cc2[C@H]([C@@H](Oc2cc1O)c1ccc(O)cc1)c1cc(O)cc(O)c1)c1ccc(O)cc1 '
              '-> False, Molecule does not contain the 1-benzofuran core\n'
              'SMILES: '
              'Oc1ccc(cc1)[C@@H]1Oc2ccc(\\C=C\\c3cc(O)cc(O)c3)cc2[C@H]1c1cc(O)cc(O)c1 '
              '-> False, Molecule does not contain the 1-benzofuran core\n'
              'SMILES: '
              'O1[C@@H]([C@H](C2=C1C(OC)=CC(=C2)C(=O)[H])C)C3=CC(OC)=C(OC)C=C3 '
              '-> False, Molecule does not contain the 1-benzofuran core\n'
              'SMILES: COc1cc(O)cc2OC(O)(Cc3ccc(O)cc3)C(=O)c12 -> False, '
              'Molecule does not contain the 1-benzofuran core\n'
              'SMILES: CC1(CC2=C(O1)C(=CC=C2)OCC(CNCC3=CC=CC=N3)O)C -> False, '
              'Molecule does not contain the 1-benzofuran core\n'
              'SMILES: '
              'C\\C=C\\c1c(O)c(O)cc2OC(=O)[C@](C)(NC(=O)\\C=C\\C(O)=O)c12 -> '
              'False, Molecule does not contain the 1-benzofuran core\n'
              'SMILES: COc1ccc2cc(oc2c1)-c1ccc(O)cc1OC -> True, Molecule is a '
              '1-benzofuran or its substituted derivative\n',
    'num_true_positives': 4,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 20,
    'precision': 1.0,
    'recall': 0.16666666666666666,
    'f1': 0.2857142857142857,
    'accuracy': None}