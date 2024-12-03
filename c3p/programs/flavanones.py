"""
Classifies: CHEBI:28863 flavanones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavanones(smiles: str):
    """
    Determines if a molecule is a flavanone (3,4-dihydro-2-aryl-2H-1-benzopyran-4-one skeleton and its substituted derivatives).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavanone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavanone core structure
    flavanone_core = "O=C1C2=C(O)C=CC=C2CC(O1)c1ccccc1"
    flavanone_core_mol = Chem.MolFromSmiles(flavanone_core)

    if flavanone_core_mol is None:
        return False, "Invalid flavanone core structure"

    # Check if the molecule contains the flavanone core
    if not mol.HasSubstructMatch(flavanone_core_mol):
        return False, "Molecule does not contain the flavanone core structure"

    return True, "Molecule contains the flavanone core structure"

# Example usage:
# smiles = "O1C(C=2C(OC)=C(C(OC)=C(O)C2O)CC=C(C)C)CC(=O)C3=C1C=C(O)C=C3"
# print(is_flavanones(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:28863',
                          'name': 'flavanones',
                          'definition': 'Members of the class of flavans with '
                                        'a '
                                        '3,4-dihydro-2-aryl-2H-1-benzopyran-4-one '
                                        'skeleton and its substituted '
                                        'derivatives.',
                          'parents': ['CHEBI:38672', 'CHEBI:3992']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 37,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}