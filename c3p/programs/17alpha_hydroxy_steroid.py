"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid (The alpha-stereoisomer of 17-hydroxy steroid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the steroid backbone
    steroid_skeleton = Chem.MolFromSmarts('C1CCC2C1CCC3C2CCC4C3(CCC(C4)O)')
    if not mol.HasSubstructMatch(steroid_skeleton):
        return False, "No steroid backbone found"

    # Check for the presence of the 17-hydroxy group with alpha stereochemistry
    hydroxyl_17alpha = Chem.MolFromSmarts('C[C@]12CC[C@@](O)(C1)CC[C@]2(C)C')
    if not mol.HasSubstructMatch(hydroxyl_17alpha):
        return False, "No 17alpha-hydroxy group found"

    return True, "Molecule is a 17alpha-hydroxy steroid"

# Example usage:
# smiles = "C[C@H]1C[C@H]2[C@@H]3CC[C@](O)(C(=O)COC(C)=O)[C@@]3(C)C[C@H](O)[C@@H]2[C@@]2(C)C=CC(=O)C=C12"
# print(is_17alpha_hydroxy_steroid(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35342',
                          'name': '17alpha-hydroxy steroid',
                          'definition': 'The alpha-stereoisomer of 17-hydroxy '
                                        'steroid.',
                          'parents': ['CHEBI:36838']},
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
    'num_true_negatives': 10,
    'num_false_negatives': 10,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}