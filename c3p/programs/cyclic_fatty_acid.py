"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a ring
    if not mol.GetRingInfo().AtomRings():
        return False, "No rings found in the molecule"

    # Check for the presence of a fatty acid group (carboxylic acid group)
    carboxylic_acid = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found in the molecule"

    # If both conditions are met, it is a cyclic fatty acid
    return True, "Molecule contains a ring and a carboxylic acid group"

# Examples for testing
test_smiles = [
    "O1C(CCCCCCCCCCCCC(O[C@H](COC(=O)CCCCCCCCCCCCC=2OC(CCCCC)=CC2C)CO)=O)=C(C(=C1CCCCC)C)C",
    "O=C(O)/C=C(/[C@H]1O[C@@H]1[C@H]2C(=C)CCCC2(C)C)\C",
    "O1C(C(O)CC)=C(C(=C1CCC(O)=O)C(O)=O)C",
    "OC(=O)CCCCCCCCCC1CCCC1",
    "O1C(=CC=C1/C=C\C(O)=O)C(=O)C#CCCCC",
    "O=C(O)/C=C/[C@@H]([C@@H]1O[C@]2(O[C@@H]([C@@H](C)[C@@H](C2)O)C)[C@@H](C)[C@@H]([C@@H]1C)O)C",
    "O=C(O)/C(=C/C[C@]12O[C@H]1[C@@H](O)C=C([C@H]2O)/C=C/CCCCC)/C",
    "C(C(O)=O)C/C=C\C/C=C\C/C=C\C/C=C\C/C=C\C[C@@H]1[C@H](CC)O1",
    "O1C(CCCCCCCCCCC(O[C@H](COC(=O)CCCCCCCCCCCCC=2OC(CCCCC)=CC2C)CO)=O)=C(C(=C1CCCCC)C)C",
    "O=C1[C@H]([C@@](C=2CC[C@]3([C@H](C2C1)CC[C@@H]3C(=O)CO)C)(CCC(=O)O)C)[C@@H](O)C",
    "O=C1N([C@@H](C[C@@](O)(C(=O)O)C(C)C)C(C1=C(O)C=CC=CCCC)=O)C",
    "C(C(O)=O)C/C=C\C/C=C\C\C=C/C=C/C=C/[C@H]1[C@H](C/C=C\CC)O1"
]

for smiles in test_smiles:
    result, reason = is_cyclic_fatty_acid(smiles)
    print(f"SMILES: {smiles} -> {result}: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:59238',
                          'name': 'cyclic fatty acid',
                          'definition': 'Any fatty acid containing anywhere in '
                                        'its structure a ring of atoms.',
                          'parents': ['CHEBI:35366']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[01:26:57] SMILES Parse Error: syntax error while parsing: '
             'O=C\x01NC(CC(O)(C(=O)O)C(C)C)C(/C1=C(\\O)/C=C/C=C/CCC)=O\n'
             '[01:26:57] SMILES Parse Error: Failed parsing SMILES '
             "'O=C\x01NC(CC(O)(C(=O)O)C(C)C)C(/C1=C(\\O)/C=C/C=C/CCC)=O' for "
             'input: '
             "'O=C\x01NC(CC(O)(C(=O)O)C(C)C)C(/C1=C(\\O)/C=C/C=C/CCC)=O'\n",
    'stdout': 'SMILES: '
              'O1C(CCCCCCCCCCCCC(O[C@H](COC(=O)CCCCCCCCCCCCC=2OC(CCCCC)=CC2C)CO)=O)=C(C(=C1CCCCC)C)C '
              '-> True: Molecule contains a ring and a carboxylic acid group\n'
              'SMILES: O=C(O)/C=C(/[C@H]1O[C@@H]1[C@H]2C(=C)CCCC2(C)C)\\C -> '
              'True: Molecule contains a ring and a carboxylic acid group\n'
              'SMILES: O1C(C(O)CC)=C(C(=C1CCC(O)=O)C(O)=O)C -> True: Molecule '
              'contains a ring and a carboxylic acid group\n'
              'SMILES: OC(=O)CCCCCCCCCC1CCCC1 -> True: Molecule contains a '
              'ring and a carboxylic acid group\n'
              'SMILES: O1C(=CC=C1/C=C\\C(O)=O)C(=O)C#CCCCC -> True: Molecule '
              'contains a ring and a carboxylic acid group\n'
              'SMILES: '
              'O=C(O)/C=C/[C@@H]([C@@H]1O[C@]2(O[C@@H]([C@@H](C)[C@@H](C2)O)C)[C@@H](C)[C@@H]([C@@H]1C)O)C '
              '-> True: Molecule contains a ring and a carboxylic acid group\n'
              'SMILES: '
              'O=C(O)/C(=C/C[C@]12O[C@H]1[C@@H](O)C=C([C@H]2O)/C=C/CCCCC)/C -> '
              'True: Molecule contains a ring and a carboxylic acid group\n'
              'SMILES: '
              'C(C(O)=O)C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C[C@@H]1[C@H](CC)O1 '
              '-> True: Molecule contains a ring and a carboxylic acid group\n'
              'SMILES: '
              'O1C(CCCCCCCCCCC(O[C@H](COC(=O)CCCCCCCCCCCCC=2OC(CCCCC)=CC2C)CO)=O)=C(C(=C1CCCCC)C)C '
              '-> True: Molecule contains a ring and a carboxylic acid group\n'
              'SMILES: '
              'O=C1[C@H]([C@@](C=2CC[C@]3([C@H](C2C1)CC[C@@H]3C(=O)CO)C)(CCC(=O)O)C)[C@@H](O)C '
              '-> True: Molecule contains a ring and a carboxylic acid group\n'
              'SMILES: '
              'O=C1N([C@@H](C[C@@](O)(C(=O)O)C(C)C)C(C1=C(O)C=CC=CCCC)=O)C -> '
              'True: Molecule contains a ring and a carboxylic acid group\n'
              'SMILES: '
              'C(C(O)=O)C/C=C\\C/C=C\\C\\C=C/C=C/C=C/[C@H]1[C@H](C/C=C\\CC)O1 '
              '-> True: Molecule contains a ring and a carboxylic acid group\n',
    'num_true_positives': 27,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9642857142857143,
    'f1': 0.9818181818181818,
    'accuracy': None}