"""
Classifies: CHEBI:27300 vitamin D
"""
from rdkit import Chem

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a member of the vitamin D class.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a vitamin D, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the secosteroid core structure
    secosteroid_core = Chem.MolFromSmarts("C1CCC2C1(CCC3C2CCC4=CC(CCC34)C)C")
    if not mol.HasSubstructMatch(secosteroid_core):
        return False, "Molecule does not contain the secosteroid core structure"

    # Check for the presence of hydroxyl groups
    hydroxyl_groups = Chem.MolFromSmarts("[OH]")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_groups))
    if hydroxyl_count < 2:
        return False, "Molecule does not have the required number of hydroxyl groups"

    return True, "Molecule is a vitamin D"

# Example usage:
# result, reason = is_vitamin_D("C1[C@]2([C@](/C(=C/C=C/3\C(CC[C@@H](C3)O)=C)/CC1)(CC[C@@]2([C@](C)(CC(C(C)C)O)O)[H])[H])C")
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:27300',
                          'name': 'vitamin D',
                          'definition': 'Any member of a group of fat-soluble '
                                        'hydroxy seco-steroids that exhibit '
                                        'biological activity against vitamin D '
                                        'deficiency. Vitamin D  can be '
                                        'obtained from sun exposure, food and '
                                        'supplements and is biologically '
                                        'inactive and converted into the '
                                        'biologically active calcitriol via '
                                        'double hydroxylation in the body.',
                          'parents': ['CHEBI:36853']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[21:16:00] SMILES Parse Error: syntax error while parsing: '
             'OC(CCC[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)CCC3=C)[H])C)[H])C)(C)C.O\n'
             '[21:16:00] SMILES Parse Error: Failed parsing SMILES '
             "'OC(CCC[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)CCC3=C)[H])C)[H])C)(C)C.O' "
             'for input: '
             "'OC(CCC[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)CCC3=C)[H])C)[H])C)(C)C.O'\n"
             '[21:16:00] SMILES Parse Error: syntax error while parsing: '
             'OC(CCC[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)[C@@H](OCCO)[C@H](O)C3=C)[H])C)[H])C)(C)C\n'
             '[21:16:00] SMILES Parse Error: Failed parsing SMILES '
             "'OC(CCC[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)[C@@H](OCCO)[C@H](O)C3=C)[H])C)[H])C)(C)C' "
             'for input: '
             "'OC(CCC[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)[C@@H](OCCO)[C@H](O)C3=C)[H])C)[H])C)(C)C'\n"
             '[21:16:00] SMILES Parse Error: syntax error while parsing: '
             'O[C@@H](\\C=C\\[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])C)C(O)(C)C\n'
             '[21:16:00] SMILES Parse Error: Failed parsing SMILES '
             "'O[C@@H](\\C=C\\[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])C)C(O)(C)C' "
             'for input: '
             "'O[C@@H](\\C=C\\[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])C)C(O)(C)C'\n"
             '[21:16:00] SMILES Parse Error: syntax error while parsing: '
             'O[C@]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])(CCCCC(O)(C)C)C\n'
             '[21:16:00] SMILES Parse Error: Failed parsing SMILES '
             "'O[C@]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])(CCCCC(O)(C)C)C' "
             'for input: '
             "'O[C@]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])(CCCCC(O)(C)C)C'\n"
             '[21:16:00] SMILES Parse Error: syntax error while parsing: '
             'OC(\\C=C\\[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)CCC3=C)[H])C)[H])C)(C(C)C)C\n'
             '[21:16:00] SMILES Parse Error: Failed parsing SMILES '
             "'OC(\\C=C\\[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)CCC3=C)[H])C)[H])C)(C(C)C)C' "
             'for input: '
             "'OC(\\C=C\\[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)CCC3=C)[H])C)[H])C)(C(C)C)C'\n"
             '[21:16:00] SMILES Parse Error: syntax error while parsing: '
             'O([C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](OC(=O)CCCC(O)=O)C[C@H](O)C3=C)[H])C)[H])C)CCC(O)(C)C\n'
             '[21:16:00] SMILES Parse Error: Failed parsing SMILES '
             "'O([C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](OC(=O)CCCC(O)=O)C[C@H](O)C3=C)[H])C)[H])C)CCC(O)(C)C' "
             'for input: '
             "'O([C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](OC(=O)CCCC(O)=O)C[C@H](O)C3=C)[H])C)[H])C)CCC(O)(C)C'\n"
             '[21:16:00] SMILES Parse Error: syntax error while parsing: '
             'FC(F)(F)C(O)(CCC[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)CCC3=C)[H])C)[H])C)C\n'
             '[21:16:00] SMILES Parse Error: Failed parsing SMILES '
             "'FC(F)(F)C(O)(CCC[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)CCC3=C)[H])C)[H])C)C' "
             'for input: '
             "'FC(F)(F)C(O)(CCC[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)CCC3=C)[H])C)[H])C)C'\n"
             '[21:16:00] SMILES Parse Error: syntax error while parsing: '
             'O[C@](\\C=C\\[C@H]([C@@H]1[C@@]2([C@@H](CC1)/C(/CCC2)=C/C=C\x03/C[C@@H](O)CCC3=C)C)C)(C(O)(C)C)C\n'
             '[21:16:00] SMILES Parse Error: Failed parsing SMILES '
             "'O[C@](\\C=C\\[C@H]([C@@H]1[C@@]2([C@@H](CC1)/C(/CCC2)=C/C=C\x03/C[C@@H](O)CCC3=C)C)C)(C(O)(C)C)C' "
             'for input: '
             "'O[C@](\\C=C\\[C@H]([C@@H]1[C@@]2([C@@H](CC1)/C(/CCC2)=C/C=C\x03/C[C@@H](O)CCC3=C)C)C)(C(O)(C)C)C'\n"
             '[21:16:00] SMILES Parse Error: syntax error while parsing: '
             'S(=O)(=O)(\\C=C\\C[C@H](C=1[C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)C)C(C)(C)C\n'
             '[21:16:00] SMILES Parse Error: Failed parsing SMILES '
             "'S(=O)(=O)(\\C=C\\C[C@H](C=1[C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)C)C(C)(C)C' "
             'for input: '
             "'S(=O)(=O)(\\C=C\\C[C@H](C=1[C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)C)C(C)(C)C'\n"
             '[21:16:00] SMILES Parse Error: syntax error while parsing: '
             'FC([C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)CCC3=C)[H])C)[H])C)CCC(C)C\n'
             '[21:16:00] SMILES Parse Error: Failed parsing SMILES '
             "'FC([C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)CCC3=C)[H])C)[H])C)CCC(C)C' "
             'for input: '
             "'FC([C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)CCC3=C)[H])C)[H])C)CCC(C)C'\n"
             '[21:16:00] SMILES Parse Error: syntax error while parsing: '
             'O[C@](\\C=C\\[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])C)(C(O)(C)C)CO\n'
             '[21:16:00] SMILES Parse Error: Failed parsing SMILES '
             "'O[C@](\\C=C\\[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])C)(C(O)(C)C)CO' "
             'for input: '
             "'O[C@](\\C=C\\[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])C)(C(O)(C)C)CO'\n"
             '[21:16:00] SMILES Parse Error: syntax error while parsing: '
             'O([C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])C)CCC(O)(CC)CC\n'
             '[21:16:00] SMILES Parse Error: Failed parsing SMILES '
             "'O([C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])C)CCC(O)(CC)CC' "
             'for input: '
             "'O([C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])C)CCC(O)(CC)CC'\n"
             '[21:16:00] SMILES Parse Error: syntax error while parsing: '
             'O[C@H]([C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)CCC3=C)[H])C)[H])C)CCC(C)C\n'
             '[21:16:00] SMILES Parse Error: Failed parsing SMILES '
             "'O[C@H]([C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)CCC3=C)[H])C)[H])C)CCC(C)C' "
             'for input: '
             "'O[C@H]([C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)CCC3=C)[H])C)[H])C)CCC(C)C'\n"
             '[21:16:00] SMILES Parse Error: syntax error while parsing: '
             'FC(F)(F)C(O)(C#CC[C@@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])CCCC(O)(C)C)C(F)(F)F\n'
             '[21:16:00] SMILES Parse Error: Failed parsing SMILES '
             "'FC(F)(F)C(O)(C#CC[C@@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])CCCC(O)(C)C)C(F)(F)F' "
             'for input: '
             "'FC(F)(F)C(O)(C#CC[C@@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\x03/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])CCCC(O)(C)C)C(F)(F)F'\n",
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 25,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}