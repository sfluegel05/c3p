"""
Classifies: CHEBI:22187 acetophenones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_acetophenones(smiles: str):
    """
    Determines if a molecule is an acetophenone (PhC(=O)CH3 and its substituted derivatives).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acetophenone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the acetophenone core structure
    acetophenone_core = Chem.MolFromSmarts("CC(=O)c1ccccc1")

    if acetophenone_core is None:
        return False, "Error in defining acetophenone core structure"

    # Check if the molecule contains the acetophenone core
    if not mol.HasSubstructMatch(acetophenone_core):
        return False, "Molecule does not contain the acetophenone core structure"

    return True, "Molecule is an acetophenone or its substituted derivative"

# Example usage
smiles_list = [
    "O1C([C@@H](O)[C@H](O)C(O)[C@@H]1OC2=C(OC)C(O)=C3C(OC=C(C3=O)C4=CC=C(O)C=C4)=C2)CO[C@@H]5OC([C@@H](O)[C@H](O)C5O)CO",
    "O=C1C(C2=CC(OC3OC(C(O)C(C3O)O)C)=C(O)C=C2)=COC4=C1C=CC(=C4)O",
    "O=C1C=C(OC)[C@@]2([C@@]1(OC=3C(=C(O)C(=C(C23)O)C)C(=O)C)O)C",
    "BrCC(=O)c1ccccc1",
    "C1(=C(C2=C(C(C(=CO2)C3=CC=C(C(=C3)OC)OC)=O)C(=C1)O[C@H]4[C@H](C([C@@H](C(O4)CO)O)O)O[C@H]5[C@H](C([C@H](C(O5)C)O)O)O)C)OC",
    "O1C(C(O)C(O)C(O)C1OC2=CC=3OC=C(C(=O)C3C(OC)=C2)C4=CC=C(O)C=C4)CO",
    "O1C([C@@H](O)[C@H](O)C(O)[C@@H]1OC2=C(OC)C=C3C(OC=C(C3=O)C4=CC(OC)=C(OC)C=C4)=C2)CO",
    "S(OC1=CC=C(C=2C(=O)C3=C(OC2)C=C(O[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)C(O)=O)C=C3)C=C1)(O)(=O)=O",
    "O=C1C(C2=CC=C(O[C@@H]3O[C@H]([C@H](O)[C@H]([C@H]3O)OC)C)C=C2)=COC4=C1C=CC(=C4)O[C@@H]5O[C@H]([C@H](O)[C@H]([C@H]5O)OC)C",
    "CC(=O)C1=CC(C)=C(O)C=C1",
    "O1C(C(O)C(O)C(O)C1OC2=C3C(OC=C(C3=O)C4=CC(OC)=C(OC)C=C4)=C(C(OC)=C2)C)CO",
    "O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1OC2=CC=3OCC(CC3C=C2)C4=C(O)C(OC)=C(OC)C=C4)CO",
    "O1C(C(O)C(O)C(O)C1OC2=CC=3OC=C(C(=O)C3C=C2)C4=CC=C(OC)C=C4)COC(=O)CC(O)=O",
    "S(OC=1C=C2OC=C(C3=CC=C(O[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)C(O)=O)C=C3)C(=O)C2=C(O)C1)(O)(=O)=O",
    "CC(C)C(=O)c1ccc(O)c(O)c1",
    "O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1OC2=CC=C(C=C2)C=3C(=O)C=4C(OC3)=CC(O)=CC4O)CO",
    "O=C1C(C2=CC=C(O)C=C2)=COC3=C1C=C(OC)C(=C3)O[C@@H]4O[C@H]([C@H](O)[C@H]([C@H]4O)O)C",
    "C=1(C(=NN(C1O)C)C)C(C2=CC=C(C=C2Cl)Cl)=O"
]

for smiles in smiles_list:
    result, reason = is_acetophenones(smiles)
    print(f"SMILES: {smiles}\nResult: {result}\nReason: {reason}\n")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22187',
                          'name': 'acetophenones',
                          'definition': 'A class or aromatic ketone consisting '
                                        'of acetophenone, PhC(=O)CH3, and its '
                                        'substituted derivatives.',
                          'parents': ['CHEBI:51867', 'CHEBI:76224']},
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
              'O1C([C@@H](O)[C@H](O)C(O)[C@@H]1OC2=C(OC)C(O)=C3C(OC=C(C3=O)C4=CC=C(O)C=C4)=C2)CO[C@@H]5OC([C@@H](O)[C@H](O)C5O)CO\n'
              'Result: False\n'
              'Reason: Molecule does not contain the acetophenone core '
              'structure\n'
              '\n'
              'SMILES: '
              'O=C1C(C2=CC(OC3OC(C(O)C(C3O)O)C)=C(O)C=C2)=COC4=C1C=CC(=C4)O\n'
              'Result: False\n'
              'Reason: Molecule does not contain the acetophenone core '
              'structure\n'
              '\n'
              'SMILES: '
              'O=C1C=C(OC)[C@@]2([C@@]1(OC=3C(=C(O)C(=C(C23)O)C)C(=O)C)O)C\n'
              'Result: True\n'
              'Reason: Molecule is an acetophenone or its substituted '
              'derivative\n'
              '\n'
              'SMILES: BrCC(=O)c1ccccc1\n'
              'Result: True\n'
              'Reason: Molecule is an acetophenone or its substituted '
              'derivative\n'
              '\n'
              'SMILES: '
              'C1(=C(C2=C(C(C(=CO2)C3=CC=C(C(=C3)OC)OC)=O)C(=C1)O[C@H]4[C@H](C([C@@H](C(O4)CO)O)O)O[C@H]5[C@H](C([C@H](C(O5)C)O)O)O)C)OC\n'
              'Result: False\n'
              'Reason: Molecule does not contain the acetophenone core '
              'structure\n'
              '\n'
              'SMILES: '
              'O1C(C(O)C(O)C(O)C1OC2=CC=3OC=C(C(=O)C3C(OC)=C2)C4=CC=C(O)C=C4)CO\n'
              'Result: False\n'
              'Reason: Molecule does not contain the acetophenone core '
              'structure\n'
              '\n'
              'SMILES: '
              'O1C([C@@H](O)[C@H](O)C(O)[C@@H]1OC2=C(OC)C=C3C(OC=C(C3=O)C4=CC(OC)=C(OC)C=C4)=C2)CO\n'
              'Result: False\n'
              'Reason: Molecule does not contain the acetophenone core '
              'structure\n'
              '\n'
              'SMILES: '
              'S(OC1=CC=C(C=2C(=O)C3=C(OC2)C=C(O[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)C(O)=O)C=C3)C=C1)(O)(=O)=O\n'
              'Result: False\n'
              'Reason: Molecule does not contain the acetophenone core '
              'structure\n'
              '\n'
              'SMILES: '
              'O=C1C(C2=CC=C(O[C@@H]3O[C@H]([C@H](O)[C@H]([C@H]3O)OC)C)C=C2)=COC4=C1C=CC(=C4)O[C@@H]5O[C@H]([C@H](O)[C@H]([C@H]5O)OC)C\n'
              'Result: False\n'
              'Reason: Molecule does not contain the acetophenone core '
              'structure\n'
              '\n'
              'SMILES: CC(=O)C1=CC(C)=C(O)C=C1\n'
              'Result: True\n'
              'Reason: Molecule is an acetophenone or its substituted '
              'derivative\n'
              '\n'
              'SMILES: '
              'O1C(C(O)C(O)C(O)C1OC2=C3C(OC=C(C3=O)C4=CC(OC)=C(OC)C=C4)=C(C(OC)=C2)C)CO\n'
              'Result: False\n'
              'Reason: Molecule does not contain the acetophenone core '
              'structure\n'
              '\n'
              'SMILES: '
              'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1OC2=CC=3OCC(CC3C=C2)C4=C(O)C(OC)=C(OC)C=C4)CO\n'
              'Result: False\n'
              'Reason: Molecule does not contain the acetophenone core '
              'structure\n'
              '\n'
              'SMILES: '
              'O1C(C(O)C(O)C(O)C1OC2=CC=3OC=C(C(=O)C3C=C2)C4=CC=C(OC)C=C4)COC(=O)CC(O)=O\n'
              'Result: False\n'
              'Reason: Molecule does not contain the acetophenone core '
              'structure\n'
              '\n'
              'SMILES: '
              'S(OC=1C=C2OC=C(C3=CC=C(O[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)C(O)=O)C=C3)C(=O)C2=C(O)C1)(O)(=O)=O\n'
              'Result: False\n'
              'Reason: Molecule does not contain the acetophenone core '
              'structure\n'
              '\n'
              'SMILES: CC(C)C(=O)c1ccc(O)c(O)c1\n'
              'Result: True\n'
              'Reason: Molecule is an acetophenone or its substituted '
              'derivative\n'
              '\n'
              'SMILES: '
              'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1OC2=CC=C(C=C2)C=3C(=O)C=4C(OC3)=CC(O)=CC4O)CO\n'
              'Result: False\n'
              'Reason: Molecule does not contain the acetophenone core '
              'structure\n'
              '\n'
              'SMILES: '
              'O=C1C(C2=CC=C(O)C=C2)=COC3=C1C=C(OC)C(=C3)O[C@@H]4O[C@H]([C@H](O)[C@H]([C@H]4O)O)C\n'
              'Result: False\n'
              'Reason: Molecule does not contain the acetophenone core '
              'structure\n'
              '\n'
              'SMILES: C=1(C(=NN(C1O)C)C)C(C2=CC=C(C=C2Cl)Cl)=O\n'
              'Result: False\n'
              'Reason: Molecule does not contain the acetophenone core '
              'structure\n'
              '\n',
    'num_true_positives': 4,
    'num_false_positives': 8,
    'num_true_negatives': 10,
    'num_false_negatives': 14,
    'precision': 0.3333333333333333,
    'recall': 0.2222222222222222,
    'f1': 0.26666666666666666,
    'accuracy': None}