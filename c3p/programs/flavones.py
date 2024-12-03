"""
Classifies: CHEBI:24043 flavones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone (2-aryl-1-benzopyran-4-one skeleton and its substituted derivatives).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavone core structure
    flavone_core = 'c1cc(-c2ccc3c(c2)oc(=O)cc3)ccc1'
    core_mol = Chem.MolFromSmiles(flavone_core)

    if core_mol is None:
        return None, "Invalid core structure"

    # Check if the molecule contains the flavone core structure
    if not mol.HasSubstructMatch(core_mol):
        return False, "Molecule does not contain the flavone core structure"

    return True, "Molecule contains the flavone core structure"

# Example test cases
smiles_list = [
    'OC[C@H]1O[C@H]([C@H](O)[C@@H](O)[C@@H]1O)c1c(O)cc(O)c2c1oc(cc2=O)-c1ccc(O)c(O)c1',
    'COc1ccc(cc1)-c1cc(=O)c2c(O)cc(O)c(-c3cc(ccc3OC)[C@@H]3CC(=O)c4c(O)cc(O)cc4O3)c2o1',
    'Oc1cc(O)c2c(c1)oc(cc2=O)-c1ccccc1',
    'COc1cc(O)c2c(c1)oc(-c1ccc(O)c(OC)c1)c(OC)c2=O'
]

for smiles in smiles_list:
    result, reason = is_flavones(smiles)
    print(f"SMILES: {smiles} => Is Flavone: {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24043',
                          'name': 'flavones',
                          'definition': 'A member of the class of flavonoid '
                                        'with a 2-aryl-1-benzopyran-4-one '
                                        '(2-arylchromen-4-one) skeleton and '
                                        'its substituted derivatives.',
                          'parents': ['CHEBI:192499', 'CHEBI:47916']},
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
              'OC[C@H]1O[C@H]([C@H](O)[C@@H](O)[C@@H]1O)c1c(O)cc(O)c2c1oc(cc2=O)-c1ccc(O)c(O)c1 '
              '=> Is Flavone: False, Reason: Molecule does not contain the '
              'flavone core structure\n'
              'SMILES: '
              'COc1ccc(cc1)-c1cc(=O)c2c(O)cc(O)c(-c3cc(ccc3OC)[C@@H]3CC(=O)c4c(O)cc(O)cc4O3)c2o1 '
              '=> Is Flavone: False, Reason: Molecule does not contain the '
              'flavone core structure\n'
              'SMILES: Oc1cc(O)c2c(c1)oc(cc2=O)-c1ccccc1 => Is Flavone: False, '
              'Reason: Molecule does not contain the flavone core structure\n'
              'SMILES: COc1cc(O)c2c(c1)oc(-c1ccc(O)c(OC)c1)c(OC)c2=O => Is '
              'Flavone: False, Reason: Molecule does not contain the flavone '
              'core structure\n',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 66,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}