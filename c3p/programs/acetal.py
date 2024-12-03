"""
Classifies: CHEBI:59769 acetal
"""
from rdkit import Chem

def is_acetal(smiles: str):
    """
    Determines if a molecule is an acetal (RR'C(OR'')(OR''') with R'', R''' != H).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acetal, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the acetal pattern
    acetal_pattern = Chem.MolFromSmarts("[#6][#6](O[#6])(O[#6])")

    if mol.HasSubstructMatch(acetal_pattern):
        return True, "Matches acetal pattern"
    else:
        return False, "Does not match acetal pattern"

# Examples
examples = [
    "O=C1C2OC(OC(C2C)C(/C=C(/C=C/C(O)=C3C(=O)NCC3=O)\C)C)(C)C4(C1)OC4C(C(=O)OC)C",
    "O=C1O[C@@H]2C[C@@]3(O[C@H](/C(=C/C)/C)[C@@H](C)[C@H](C3)O)O[C@@H](C2)CC=C(C[C@H](C=CC=C4[C@]5([C@H]1C=C([C@@H](OC)[C@H]5OC4)C)O)C)C",
    "CCCN(CC)CC1COC2(CCC(CC2)C(C)(C)C)O1",
    "O1C2(OCCCC2)CCC1C"
]

for example in examples:
    result, reason = is_acetal(example)
    print(f"SMILES: {example} -> {result}, {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:59769',
                          'name': 'acetal',
                          'definition': 'An organooxygen compound having the '
                                        "structure RR'C(OR'')(OR''') (R'', "
                                        "R''' =/= H). Mixed acetals have R'' "
                                        "and R''' groups which differ.",
                          'parents': ['CHEBI:36963']},
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
              'O=C1C2OC(OC(C2C)C(/C=C(/C=C/C(O)=C3C(=O)NCC3=O)\\C)C)(C)C4(C1)OC4C(C(=O)OC)C '
              '-> True, Matches acetal pattern\n'
              'SMILES: '
              'O=C1O[C@@H]2C[C@@]3(O[C@H](/C(=C/C)/C)[C@@H](C)[C@H](C3)O)O[C@@H](C2)CC=C(C[C@H](C=CC=C4[C@]5([C@H]1C=C([C@@H](OC)[C@H]5OC4)C)O)C)C '
              '-> True, Matches acetal pattern\n'
              'SMILES: CCCN(CC)CC1COC2(CCC(CC2)C(C)(C)C)O1 -> True, Matches '
              'acetal pattern\n'
              'SMILES: O1C2(OCCCC2)CCC1C -> True, Matches acetal pattern\n',
    'num_true_positives': 41,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 7,
    'precision': 0.9761904761904762,
    'recall': 0.8541666666666666,
    'f1': 0.9111111111111111,
    'accuracy': None}