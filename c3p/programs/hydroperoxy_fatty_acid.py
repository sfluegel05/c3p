"""
Classifies: CHEBI:64009 hydroperoxy fatty acid
"""
from rdkit import Chem

def is_hydroperoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroperoxy fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroperoxy fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Check for the presence of hydroperoxy groups
    hydroperoxy = Chem.MolFromSmarts('OO')
    if not mol.HasSubstructMatch(hydroperoxy):
        return False, "No hydroperoxy group found"

    # Check for long carbon chain (fatty acid characteristic)
    carbon_chain = Chem.MolFromSmarts('[C;R0]-[C;R0]-[C;R0]-[C;R0]-[C;R0]')
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "No long carbon chain found"

    return True, "Molecule is a hydroperoxy fatty acid"

# Example usage
smiles_examples = [
    "OC(CCCCC/C=C\C/C=C\C=C\C(C/C=C\C/C=C\CC)OO)=O",
    "O1OC\2CC1C(/C2=C/C(OO)CCCC(O)=O)C/C=C/CCCCC",
    "CCCCC\C=C/C[C@@H](OO)\C=C\CCCCCCC(O)=O",
    "CC\C=C/C\C=C/C[C@H](OO)\C=C\C=C/C\C=C/C\C=C/CCC(O)=O",
    "[C@@H]1([C@@H]2C[C@H]([C@@H]1/C=C/[C@H](C/C=C\CC)OO)OO2)C/C=C\CCCC(O)=O",
    "C(CCC)CC(/C=C/C=C\C\C=C/C=C/C(CCCC(O)=O)OO)O",
    "OC(CC/C=C\C/C=C\C/C=C\C/C=C\C=C\[C@H](CCCCC)OO)=O",
    "CCCCC\C=C/C[C@H](OO)\C=C\C=C/C\C=C/CCCC(O)=O",
    "C(CCC(O)=O)/C=C\C[C@@H](/C=C/C=C\C[C@H]1[C@@H](CCCCC)O1)OO",
    "C(C(O)=O)(*)OO",
    "C(CCC(O)=O)/C=C\C/C=C\C=C\[C@H](C[C@H]1[C@@H](CCCCC)O1)OO",
    "C(C[C@H]1[C@@H]2C[C@H]([C@@H]1/C=C/[C@H](CCCCC)OO)OO2)CCCCC(O)=O"
]

for smiles in smiles_examples:
    result, reason = is_hydroperoxy_fatty_acid(smiles)
    print(f"SMILES: {smiles} -> {result}, {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:64009',
                          'name': 'hydroperoxy fatty acid',
                          'definition': 'Any fatty acid carrying one or more '
                                        'hydroperoxy substituents.',
                          'parents': ['CHEBI:35366', 'CHEBI:61051']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[02:15:18] SMILES Parse Error: syntax error while parsing: '
             'O1OC\x02CC1C(/C2=C/C(OO)CCCC(O)=O)C/C=C/CCCCC\n'
             '[02:15:18] SMILES Parse Error: Failed parsing SMILES '
             "'O1OC\x02CC1C(/C2=C/C(OO)CCCC(O)=O)C/C=C/CCCCC' for input: "
             "'O1OC\x02CC1C(/C2=C/C(OO)CCCC(O)=O)C/C=C/CCCCC'\n"
             '[02:15:18] SMILES Parse Error: syntax error while parsing: '
             'O1OC\x02CC1C(/C2=C/C(OO)CCCC(O)=O)C/C=C/CCCCC\n'
             '[02:15:18] SMILES Parse Error: Failed parsing SMILES '
             "'O1OC\x02CC1C(/C2=C/C(OO)CCCC(O)=O)C/C=C/CCCCC' for input: "
             "'O1OC\x02CC1C(/C2=C/C(OO)CCCC(O)=O)C/C=C/CCCCC'\n",
    'stdout': 'SMILES: OC(CCCCC/C=C\\C/C=C\\C=C\\C(C/C=C\\C/C=C\\CC)OO)=O -> '
              'True, Molecule is a hydroperoxy fatty acid\n'
              'SMILES: O1OC\x02CC1C(/C2=C/C(OO)CCCC(O)=O)C/C=C/CCCCC -> False, '
              'Invalid SMILES string\n'
              'SMILES: CCCCC\\C=C/C[C@@H](OO)\\C=C\\CCCCCCC(O)=O -> True, '
              'Molecule is a hydroperoxy fatty acid\n'
              'SMILES: '
              'CC\\C=C/C\\C=C/C[C@H](OO)\\C=C\\C=C/C\\C=C/C\\C=C/CCC(O)=O -> '
              'False, No long carbon chain found\n'
              'SMILES: '
              '[C@@H]1([C@@H]2C[C@H]([C@@H]1/C=C/[C@H](C/C=C\\CC)OO)OO2)C/C=C\\CCCC(O)=O '
              '-> True, Molecule is a hydroperoxy fatty acid\n'
              'SMILES: C(CCC)CC(/C=C/C=C\\C\\C=C/C=C/C(CCCC(O)=O)OO)O -> True, '
              'Molecule is a hydroperoxy fatty acid\n'
              'SMILES: OC(CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C=C\\[C@H](CCCCC)OO)=O '
              '-> True, Molecule is a hydroperoxy fatty acid\n'
              'SMILES: CCCCC\\C=C/C[C@H](OO)\\C=C\\C=C/C\\C=C/CCCC(O)=O -> '
              'True, Molecule is a hydroperoxy fatty acid\n'
              'SMILES: '
              'C(CCC(O)=O)/C=C\\C[C@@H](/C=C/C=C\\C[C@H]1[C@@H](CCCCC)O1)OO -> '
              'True, Molecule is a hydroperoxy fatty acid\n'
              'SMILES: C(C(O)=O)(*)OO -> False, No long carbon chain found\n'
              'SMILES: '
              'C(CCC(O)=O)/C=C\\C/C=C\\C=C\\[C@H](C[C@H]1[C@@H](CCCCC)O1)OO -> '
              'True, Molecule is a hydroperoxy fatty acid\n'
              'SMILES: '
              'C(C[C@H]1[C@@H]2C[C@H]([C@@H]1/C=C/[C@H](CCCCC)OO)OO2)CCCCC(O)=O '
              '-> True, Molecule is a hydroperoxy fatty acid\n',
    'num_true_positives': 9,
    'num_false_positives': 0,
    'num_true_negatives': 12,
    'num_false_negatives': 3,
    'precision': 1.0,
    'recall': 0.75,
    'f1': 0.8571428571428571,
    'accuracy': None}