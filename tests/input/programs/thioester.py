"""
Classifies: CHEBI:51277 thioester
"""
from rdkit import Chem

def is_thioester(smiles: str):
    """
    Determines if a molecule is a thioester (RC(=O)SR).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thioester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define thioester SMARTS pattern
    thioester_smarts = "[#6][C](=[O])[S][#6]"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)

    # Check if the molecule matches the thioester pattern
    if mol.HasSubstructMatch(thioester_pattern):
        return True, "Molecule is a thioester"
    else:
        return False, "Molecule does not match thioester pattern"

# Example usage
smiles_list = [
    "CC(=CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)C(C)(C)C=CC(O)=O",
    "CC(C)CCCCCCCCCCCC(O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12",
    "O1C(CCC1)C(=O)N2CCNCC2"
]

for smiles in smiles_list:
    result, reason = is_thioester(smiles)
    print(f"SMILES: {smiles}, Is Thioester: {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51277',
                          'name': 'thioester',
                          'definition': 'A compound of general formula '
                                        "RC(=O)SR'. Compare with thionoester, "
                                        "RC(=S)OR'.",
                          'parents': ['CHEBI:26959', 'CHEBI:36586']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[01:07:11] SMILES Parse Error: syntax error while parsing: '
             'O=C(N[C@H]1[C@H](O[C@H](CO)[C@H]([C@@H]1O)O)OC/C(=C\x02/C3C(C(C)CC3)CC(C2)(C)C)/C)C\n'
             '[01:07:11] SMILES Parse Error: Failed parsing SMILES '
             "'O=C(N[C@H]1[C@H](O[C@H](CO)[C@H]([C@@H]1O)O)OC/C(=C\x02/C3C(C(C)CC3)CC(C2)(C)C)/C)C' "
             'for input: '
             "'O=C(N[C@H]1[C@H](O[C@H](CO)[C@H]([C@@H]1O)O)OC/C(=C\x02/C3C(C(C)CC3)CC(C2)(C)C)/C)C'\n",
    'stdout': 'SMILES: '
              'CC(=CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)C(C)(C)C=CC(O)=O, '
              'Is Thioester: True, Reason: Molecule is a thioester\n'
              'SMILES: '
              'CC(C)CCCCCCCCCCCC(O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12, '
              'Is Thioester: True, Reason: Molecule is a thioester\n'
              'SMILES: O1C(CCC1)C(=O)N2CCNCC2, Is Thioester: False, Reason: '
              'Molecule does not match thioester pattern\n',
    'num_true_positives': 82,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 3,
    'precision': 0.9879518072289156,
    'recall': 0.9647058823529412,
    'f1': 0.9761904761904762,
    'accuracy': None}