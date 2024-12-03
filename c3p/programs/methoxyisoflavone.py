"""
Classifies: CHEBI:38756 methoxyisoflavone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_methoxyisoflavone(smiles: str):
    """
    Determines if a molecule is a methoxyisoflavone (isoflavone with at least one methoxy substituent).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methoxyisoflavone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the isoflavone core structure
    isoflavone_smarts = 'O=C1C=CC2=CC=CC=C2C1C3=CC=CC=C3'
    isoflavone_core = Chem.MolFromSmarts(isoflavone_smarts)
    
    if not mol.HasSubstructMatch(isoflavone_core):
        return False, "No isoflavone core structure found"

    # Check for methoxy groups (OCH3)
    methoxy_smarts = 'CO'
    methoxy_group = Chem.MolFromSmarts(methoxy_smarts)
    
    if not mol.HasSubstructMatch(methoxy_group):
        return False, "No methoxy groups found"

    return True, "Methoxyisoflavone identified"

# Examples to test the function
smiles_examples = [
    "O1C2=C(C(=O)C(C3=C(O)C=C(OC)C=C3)=C1)C(O)=CC(O)=C2OC",  # 5,7-Dihydroxy-8,4'-dimethoxyisoflavone
    "O1C=2C(C(=O)C(C=3C(O)=C(CC=C(C)C)C(OC)=C(OC)C3)=C1)=C(O)C=C(O)C2",  # 2'-Hydroxypiscerythrinetin
    "O1C2=C(CC=C(C)C)C(O)=CC(O)=C2C(=O)C(C3=CC=C(OC)C=C3)=C1",  # Gancaonin M
    "COc1cc(O)c(c(OC)c1CC=C(C)C)-c1coc2cc(O)ccc2c1=O",  # licoricone
    "COc1c(N)c(O)c(CC=C(C)C)c(c1CC=C(C)C)-c1coc2cc(O)cc(O)c2c1=O",  # piscerythramine
    "COC1=CC=C(C=C1)C2=COC3=C(C2=O)C=CC(=C3OC)O",  # 7-hydroxy-8-methoxy-3-(4-methoxyphenyl)-1-benzopyran-4-one
    "O1C2=C(C(=O)C(C=3C(O)=CC(OC)=C(OC)C3)=C1)C=CC(O)=C2",  # 2',7-Dihydroxy-4',5'-dimethoxyisoflavone
    "COc1ccc2c(c1)occ(-c1ccc(O)c(OC)c1)c2=O",  # Sayanedine
    "O1C2=C(C(=O)C(C3=C(OC)C=C(OC)C=C3)=C1)C=CC(O)=C2",  # 2'-Methoxyformonetin
    "O1C=2C(C(=O)C(C3=CC=C(OC)C=C3)=C1)=CC(O)=C(OC)C2"  # Alfalone
]

for smiles in smiles_examples:
    result, reason = is_methoxyisoflavone(smiles)
    print(f"SMILES: {smiles} -> {result}, {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38756',
                          'name': 'methoxyisoflavone',
                          'definition': 'Members of the class of isoflavones '
                                        'with at least one methoxy '
                                        'substituent.',
                          'parents': ['CHEBI:25698', 'CHEBI:38757']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: O1C2=C(C(=O)C(C3=C(O)C=C(OC)C=C3)=C1)C(O)=CC(O)=C2OC -> '
              'False, No isoflavone core structure found\n'
              'SMILES: '
              'O1C=2C(C(=O)C(C=3C(O)=C(CC=C(C)C)C(OC)=C(OC)C3)=C1)=C(O)C=C(O)C2 '
              '-> False, No isoflavone core structure found\n'
              'SMILES: O1C2=C(CC=C(C)C)C(O)=CC(O)=C2C(=O)C(C3=CC=C(OC)C=C3)=C1 '
              '-> False, No isoflavone core structure found\n'
              'SMILES: COc1cc(O)c(c(OC)c1CC=C(C)C)-c1coc2cc(O)ccc2c1=O -> '
              'False, No isoflavone core structure found\n'
              'SMILES: '
              'COc1c(N)c(O)c(CC=C(C)C)c(c1CC=C(C)C)-c1coc2cc(O)cc(O)c2c1=O -> '
              'False, No isoflavone core structure found\n'
              'SMILES: COC1=CC=C(C=C1)C2=COC3=C(C2=O)C=CC(=C3OC)O -> False, No '
              'isoflavone core structure found\n'
              'SMILES: O1C2=C(C(=O)C(C=3C(O)=CC(OC)=C(OC)C3)=C1)C=CC(O)=C2 -> '
              'False, No isoflavone core structure found\n'
              'SMILES: COc1ccc2c(c1)occ(-c1ccc(O)c(OC)c1)c2=O -> False, No '
              'isoflavone core structure found\n'
              'SMILES: O1C2=C(C(=O)C(C3=C(OC)C=C(OC)C=C3)=C1)C=CC(O)=C2 -> '
              'False, No isoflavone core structure found\n'
              'SMILES: O1C=2C(C(=O)C(C3=CC=C(OC)C=C3)=C1)=CC(O)=C(OC)C2 -> '
              'False, No isoflavone core structure found\n',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 10,
    'num_false_negatives': 10,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}