"""
Classifies: CHEBI:26377 pterocarpans
"""
from rdkit import Chem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pterocarpan, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the pterocarpan core structure
    pterocarpan_smarts = 'C1Oc2ccccc2C2Oc3ccccc3C12'
    pterocarpan_pattern = Chem.MolFromSmarts(pterocarpan_smarts)

    if pterocarpan_pattern is None:
        return False, "Invalid SMARTS pattern for pterocarpan"

    if mol.HasSubstructMatch(pterocarpan_pattern):
        return True, "Molecule matches the pterocarpan core structure"
    else:
        return False, "Molecule does not match the pterocarpan core structure"

# Example usage
smiles_examples = [
    "O1[C@@]2([C@](O)(C3=C1C(C[C@H](O)C(O)(C)C)=C(OC)C=C3)COC4=C2C=CC(O)=C4)[H]",
    "[H][C@]12COc3cc(O)ccc3[C@@]1([H])Oc1cc(OC)ccc21",
    "O1C2C(C3=C1C4=C(OC(C=C4)(C)C)C=C3)COC5=C2C=CC(O)=C5OC",
    "[H][C@@]12Oc3cc(O)ccc3[C@]1(O)COc1cc(O)c(CC=C(C)C)cc21",
    "O1C2C(C=3C1=CC(OC)=C(O)C3)COC4=C2C=CC(OC)=C4",
    "COc1c2O[C@@H]3[C@@H](COc4cc(O)ccc34)c2ccc1O",
    "OC[C@H]1O[C@@H](Oc2ccc3[C@@H]4Oc5cc6OCOc6cc5[C@@H]4COc3c2)[C@H](O)[C@@H](O)[C@@H]1O",
    "COc1ccc2c(O[C@H]3c4ccc(O)cc4OC[C@@]23O)c1CC=C(C)C",
    "O1C2C(C=3C1=C(OC)C(OC)=C(O)C3)COC=4C2=CC(O)=C(OC)C4OC",
    "C1Oc2ccccc2C2Oc3ccccc3C12"
]

for smiles in smiles_examples:
    result, reason = is_pterocarpans(smiles)
    print(f"SMILES: {smiles} -> {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26377',
                          'name': 'pterocarpans',
                          'definition': 'Members of the class of '
                                        'benzofurochromene with a '
                                        '6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene '
                                        'skeleton and its substituted '
                                        'derivatives. They generally bear '
                                        'structural resemblance to '
                                        'isoflavanoids that possess antibiotic '
                                        'activity and are produced by plant '
                                        'tissues in response to infection. '
                                        'They are the 3,4-dihydroderivatives '
                                        'of coumestans.',
                          'parents': ['CHEBI:38834', 'CHEBI:72544']},
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
              'O1[C@@]2([C@](O)(C3=C1C(C[C@H](O)C(O)(C)C)=C(OC)C=C3)COC4=C2C=CC(O)=C4)[H] '
              '-> True, Reason: Molecule matches the pterocarpan core '
              'structure\n'
              'SMILES: [H][C@]12COc3cc(O)ccc3[C@@]1([H])Oc1cc(OC)ccc21 -> '
              'True, Reason: Molecule matches the pterocarpan core structure\n'
              'SMILES: O1C2C(C3=C1C4=C(OC(C=C4)(C)C)C=C3)COC5=C2C=CC(O)=C5OC '
              '-> True, Reason: Molecule matches the pterocarpan core '
              'structure\n'
              'SMILES: [H][C@@]12Oc3cc(O)ccc3[C@]1(O)COc1cc(O)c(CC=C(C)C)cc21 '
              '-> True, Reason: Molecule matches the pterocarpan core '
              'structure\n'
              'SMILES: O1C2C(C=3C1=CC(OC)=C(O)C3)COC4=C2C=CC(OC)=C4 -> True, '
              'Reason: Molecule matches the pterocarpan core structure\n'
              'SMILES: COc1c2O[C@@H]3[C@@H](COc4cc(O)ccc34)c2ccc1O -> True, '
              'Reason: Molecule matches the pterocarpan core structure\n'
              'SMILES: '
              'OC[C@H]1O[C@@H](Oc2ccc3[C@@H]4Oc5cc6OCOc6cc5[C@@H]4COc3c2)[C@H](O)[C@@H](O)[C@@H]1O '
              '-> True, Reason: Molecule matches the pterocarpan core '
              'structure\n'
              'SMILES: COc1ccc2c(O[C@H]3c4ccc(O)cc4OC[C@@]23O)c1CC=C(C)C -> '
              'True, Reason: Molecule matches the pterocarpan core structure\n'
              'SMILES: O1C2C(C=3C1=C(OC)C(OC)=C(O)C3)COC=4C2=CC(O)=C(OC)C4OC '
              '-> True, Reason: Molecule matches the pterocarpan core '
              'structure\n'
              'SMILES: C1Oc2ccccc2C2Oc3ccccc3C12 -> True, Reason: Molecule '
              'matches the pterocarpan core structure\n',
    'num_true_positives': 10,
    'num_false_positives': 0,
    'num_true_negatives': 10,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}