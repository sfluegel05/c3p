"""
Classifies: CHEBI:72572 isoflavans
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_isoflavans(smiles: str):
    """
    Determines if a molecule is an isoflavan (Any isoflavonoid with a 3,4-dihydro-3-aryl-2H-1-benzopyran skeleton and its substituted derivatives).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavan, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the 3,4-dihydro-2H-1-benzopyran skeleton
    benzopyran_pattern = Chem.MolFromSmarts('C1CC2=CC=CC=C2OC1')
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No 3,4-dihydro-2H-1-benzopyran skeleton found"

    # Check for the 3-aryl substitution
    aryl_substitution_pattern = Chem.MolFromSmarts('C1CC2=CC=CC=C2OC1-c3ccccc3')
    if not mol.HasSubstructMatch(aryl_substitution_pattern):
        return False, "No 3-aryl substitution found"

    return True, "Molecule is an isoflavan"

# Example usage
smiles_list = [
    "O1[C@@]2([C@H](O)C[C@H]3C(CCC[C@@]3([C@H]2C[C@@]4([C@H]1[C@H](O)C(=C)[C@@H]([C@H]4O)O)O)C)(C)C)C",
    "O=C1C2=C(OC[C@@]1(O)C3=CC=C(O)C=C3)C=C(O)C=C2",
    "CC(=C)C1Cc2c(O1)cc1oc3oc4cc(O)ccc4c3c(=O)c1c2O",
    "OC1Oc2cc(O)cc(O)c2C(=O)C1c1ccc(O)cc1",
    "O1CC([C@@H](O)C=2C1=CC=CC2)C=3C(O)=CC=CC3",
    "O1CC(C(=O)C2=C1C(O)=C(O)C=C2)C3=CC=C(O)C=C3",
    "O1C(C=CC=2C1=C(CC(O)C(C)=C)C=3OC=C(C(=O)C3C2O)C4=CC=C(O)C=C4)(C)C",
    "COc1ccc([C@H]2COc3cc(O)ccc3C2)c(O)c1",
    "C12=CC=C(C=C2OC[C@H](C1O)C=3C(=CC(=CC3)OC)O)O",
    "CC(C)=CCc1c(O)cc(O)c2C(=O)C(COc12)c1ccc(O)cc1O",
    "O1CC(C2=C(O)C(=C(O)C=C2)CC=C(C)C)C(=O)C=3C1=CC(OC)=C(C3O)CC=C(C)C",
    "O=C1C2=C(OC[C@H]1C3=CC(O)=C(O)C=C3)C(=C(O)C=C2O)C",
    "CC1(C)Oc2ccc([C@@H]3COc4cc(O)ccc4C3)c(O)c2C=C1",
    "S(OC1=C(C=C(C2COC=3C(C2=O)=C(O)C=C(O)C3CC=C(C)C)C=C1O)CC=C(C)C)(O)(=O)=O",
    "COc1ccc([C@@H]2COc3cc(O)ccc3C2=O)c(O)c1",
    "O1C(CCC=C(C)C)(C=CC=2C1=CC=3OCC(C(=O)C3C2O)C4=C(O)C=C(O)C=C4)C",
    "O1C(C(OCCO)(C)C)CC=2C1=C(C=3OC=C(C(=O)C3C2O)C4=CC=C(O)C=C4)CC=C(C)C"
]

for smiles in smiles_list:
    result, reason = is_isoflavans(smiles)
    print(f"SMILES: {smiles}\nResult: {result}, Reason: {reason}\n")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:72572',
                          'name': 'isoflavans',
                          'definition': 'Any isoflavonoid with a '
                                        '3,4-dihydro-3-aryl-2H-1-benzopyran '
                                        'skeleton and its substituted '
                                        'derivatives.',
                          'parents': ['CHEBI:50753']},
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
              'O1[C@@]2([C@H](O)C[C@H]3C(CCC[C@@]3([C@H]2C[C@@]4([C@H]1[C@H](O)C(=C)[C@@H]([C@H]4O)O)O)C)(C)C)C\n'
              'Result: False, Reason: No 3,4-dihydro-2H-1-benzopyran skeleton '
              'found\n'
              '\n'
              'SMILES: O=C1C2=C(OC[C@@]1(O)C3=CC=C(O)C=C3)C=C(O)C=C2\n'
              'Result: False, Reason: No 3,4-dihydro-2H-1-benzopyran skeleton '
              'found\n'
              '\n'
              'SMILES: CC(=C)C1Cc2c(O1)cc1oc3oc4cc(O)ccc4c3c(=O)c1c2O\n'
              'Result: False, Reason: No 3,4-dihydro-2H-1-benzopyran skeleton '
              'found\n'
              '\n'
              'SMILES: OC1Oc2cc(O)cc(O)c2C(=O)C1c1ccc(O)cc1\n'
              'Result: False, Reason: No 3,4-dihydro-2H-1-benzopyran skeleton '
              'found\n'
              '\n'
              'SMILES: O1CC([C@@H](O)C=2C1=CC=CC2)C=3C(O)=CC=CC3\n'
              'Result: False, Reason: No 3,4-dihydro-2H-1-benzopyran skeleton '
              'found\n'
              '\n'
              'SMILES: O1CC(C(=O)C2=C1C(O)=C(O)C=C2)C3=CC=C(O)C=C3\n'
              'Result: False, Reason: No 3,4-dihydro-2H-1-benzopyran skeleton '
              'found\n'
              '\n'
              'SMILES: '
              'O1C(C=CC=2C1=C(CC(O)C(C)=C)C=3OC=C(C(=O)C3C2O)C4=CC=C(O)C=C4)(C)C\n'
              'Result: False, Reason: No 3,4-dihydro-2H-1-benzopyran skeleton '
              'found\n'
              '\n'
              'SMILES: COc1ccc([C@H]2COc3cc(O)ccc3C2)c(O)c1\n'
              'Result: False, Reason: No 3,4-dihydro-2H-1-benzopyran skeleton '
              'found\n'
              '\n'
              'SMILES: C12=CC=C(C=C2OC[C@H](C1O)C=3C(=CC(=CC3)OC)O)O\n'
              'Result: False, Reason: No 3,4-dihydro-2H-1-benzopyran skeleton '
              'found\n'
              '\n'
              'SMILES: CC(C)=CCc1c(O)cc(O)c2C(=O)C(COc12)c1ccc(O)cc1O\n'
              'Result: False, Reason: No 3,4-dihydro-2H-1-benzopyran skeleton '
              'found\n'
              '\n'
              'SMILES: '
              'O1CC(C2=C(O)C(=C(O)C=C2)CC=C(C)C)C(=O)C=3C1=CC(OC)=C(C3O)CC=C(C)C\n'
              'Result: False, Reason: No 3,4-dihydro-2H-1-benzopyran skeleton '
              'found\n'
              '\n'
              'SMILES: O=C1C2=C(OC[C@H]1C3=CC(O)=C(O)C=C3)C(=C(O)C=C2O)C\n'
              'Result: False, Reason: No 3,4-dihydro-2H-1-benzopyran skeleton '
              'found\n'
              '\n'
              'SMILES: CC1(C)Oc2ccc([C@@H]3COc4cc(O)ccc4C3)c(O)c2C=C1\n'
              'Result: False, Reason: No 3,4-dihydro-2H-1-benzopyran skeleton '
              'found\n'
              '\n'
              'SMILES: '
              'S(OC1=C(C=C(C2COC=3C(C2=O)=C(O)C=C(O)C3CC=C(C)C)C=C1O)CC=C(C)C)(O)(=O)=O\n'
              'Result: False, Reason: No 3,4-dihydro-2H-1-benzopyran skeleton '
              'found\n'
              '\n'
              'SMILES: COc1ccc([C@@H]2COc3cc(O)ccc3C2=O)c(O)c1\n'
              'Result: False, Reason: No 3,4-dihydro-2H-1-benzopyran skeleton '
              'found\n'
              '\n'
              'SMILES: '
              'O1C(CCC=C(C)C)(C=CC=2C1=CC=3OCC(C(=O)C3C2O)C4=C(O)C=C(O)C=C4)C\n'
              'Result: False, Reason: No 3,4-dihydro-2H-1-benzopyran skeleton '
              'found\n'
              '\n'
              'SMILES: '
              'O1C(C(OCCO)(C)C)CC=2C1=C(C=3OC=C(C(=O)C3C2O)C4=CC=C(O)C=C4)CC=C(C)C\n'
              'Result: False, Reason: No 3,4-dihydro-2H-1-benzopyran skeleton '
              'found\n'
              '\n',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 17,
    'num_false_negatives': 17,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}