"""
Classifies: CHEBI:79202 ascaroside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_ascaroside(smiles: str):
    """
    Determines if a molecule is an ascaroside (Any glycoside derived from alpha-ascarylopyranose or its derivatives).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ascaroside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of alpha-ascarylopyranose core
    alpha_ascarylopyranose_smiles = "C[C@H]1O[C@@H](O)[C@@H](C)[C@H](O)[C@H]1O"
    alpha_ascarylopyranose = Chem.MolFromSmiles(alpha_ascarylopyranose_smiles)
    if alpha_ascarylopyranose is None:
        return False, "Error in defining alpha-ascarylopyranose structure"

    # Check if the molecule contains the alpha-ascarylopyranose core
    if not mol.HasSubstructMatch(alpha_ascarylopyranose):
        return False, "Molecule does not contain alpha-ascarylopyranose core"

    # Check for glycosidic linkage and long-chain fatty acid or derivative
    glycosidic_linkage_smiles = "O[C@@H]1O[C@@H](C)[C@H](O)[C@H](O)[C@H]1O"
    glycosidic_linkage = Chem.MolFromSmiles(glycosidic_linkage_smiles)
    if glycosidic_linkage is None:
        return False, "Error in defining glycosidic linkage structure"

    if not mol.HasSubstructMatch(glycosidic_linkage):
        return False, "Molecule does not contain glycosidic linkage"

    # Check for the presence of a long-chain fatty acid or derivative
    fatty_acid_smiles = "CCCCCCCCCCCCC(O)=O"
    fatty_acid = Chem.MolFromSmiles(fatty_acid_smiles)
    if fatty_acid is None:
        return False, "Error in defining long-chain fatty acid structure"

    if not mol.HasSubstructMatch(fatty_acid):
        return False, "Molecule does not contain long-chain fatty acid or derivative"

    return True, "Molecule is an ascaroside"

# Example usage:
# smiles = "C[C@H](CCCCCCCCCCCCC(O)=O)O[C@@H]1O[C@@H](C)[C@H](O)C[C@H]1O"
# print(is_ascaroside(smiles))  # Expected output: (True, "Molecule is an ascaroside")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:79202',
                          'name': 'ascaroside',
                          'definition': 'Any glycoside derived from '
                                        'alpha-ascarylopyranose or its '
                                        'derivatives.',
                          'parents': ['CHEBI:24400', 'CHEBI:35315']},
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
    'num_true_negatives': 15,
    'num_false_negatives': 15,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}