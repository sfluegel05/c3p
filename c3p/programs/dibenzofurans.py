"""
Classifies: CHEBI:38922 dibenzofurans
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_dibenzofurans(smiles: str):
    """
    Determines if a molecule is a dibenzofuran or its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dibenzofuran, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for dibenzofuran
    dibenzofuran_smarts = "c1ccc2c(c1)oc3c2cccc3"
    dibenzofuran_pattern = Chem.MolFromSmarts(dibenzofuran_smarts)
    
    if mol.HasSubstructMatch(dibenzofuran_pattern):
        return True, "Molecule is a dibenzofuran or its substituted derivative"
    else:
        return False, "Molecule does not match dibenzofuran structure"

# Test cases
print(is_dibenzofurans("O=C(O)C1=C(O)C=C2OC=3C(C2=C1CCCCC)=C(C=C(OC)C3)CCCCC"))
print(is_dibenzofurans("ClC1=C(O)C=C(C)C2=C1OC=3C(Cl)=C(O)C=C(C23)C"))
print(is_dibenzofurans("Oc1ccc2oc3c(Cl)c(Cl)c(Cl)cc3c2c1"))
print(is_dibenzofurans("Oc1ccc(Cl)c2oc3ccc(Cl)cc3c12"))
print(is_dibenzofurans("O1C2=C(C3=C(O)C=C4OC=5C(C4=C3C)=C(C=C(O)C5)C)C(O)=CC(=C2C6=C1C=C(O)C=C6C)C"))
print(is_dibenzofurans("Clc1cc2oc3cc(Cl)c(Cl)c(Cl)c3c2c(Cl)c1Cl"))
print(is_dibenzofurans("COC1=C(C=C2C(=C1)C3=CC=CC=C3O2)NC(=O)NC4=CC=C(C=C4)F"))
print(is_dibenzofurans("COc1cc2oc3c(OC)c(OC)c(O)cc3c2cc1O"))
print(is_dibenzofurans("Oc1cccc2oc3c(Cl)c(Cl)ccc3c12"))
print(is_dibenzofurans("COc1c(cc(O)c2c1oc1cc(O)c(O)cc21)-c1ccc(O)c(O)c1"))
print(is_dibenzofurans("O1C2=C(C(=CC(=C2)O)CO)C3=C1C=C(O)C=C3C"))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38922',
                          'name': 'dibenzofurans',
                          'definition': 'Any organic heterotricyclic compound '
                                        'based on a dibenzofuran skeleton and '
                                        'its substituted derivatives thereof.',
                          'parents': ['CHEBI:26979', 'CHEBI:38104']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(True, 'Molecule is a dibenzofuran or its substituted "
              "derivative')\n"
              "(True, 'Molecule is a dibenzofuran or its substituted "
              "derivative')\n"
              "(True, 'Molecule is a dibenzofuran or its substituted "
              "derivative')\n"
              "(True, 'Molecule is a dibenzofuran or its substituted "
              "derivative')\n"
              "(True, 'Molecule is a dibenzofuran or its substituted "
              "derivative')\n"
              "(True, 'Molecule is a dibenzofuran or its substituted "
              "derivative')\n"
              "(True, 'Molecule is a dibenzofuran or its substituted "
              "derivative')\n"
              "(True, 'Molecule is a dibenzofuran or its substituted "
              "derivative')\n"
              "(True, 'Molecule is a dibenzofuran or its substituted "
              "derivative')\n"
              "(True, 'Molecule is a dibenzofuran or its substituted "
              "derivative')\n"
              "(True, 'Molecule is a dibenzofuran or its substituted "
              "derivative')\n",
    'num_true_positives': 11,
    'num_false_positives': 0,
    'num_true_negatives': 11,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}