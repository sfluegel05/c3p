"""
Classifies: CHEBI:50753 isoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is a 1-benzopyran with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Define a flexible 1-benzopyran core
    # [O;X2] is oxygen with exactly two bonds,
    # [C;X4,X3] is a carbon with either 3 or 4 bonds, allowing for C=O or CH2 in position 4
    # c1ccccc2 is the benzene ring in the benzopyran system
    benzopyran_pattern = Chem.MolFromSmarts("c1ccccc2[O;X2][C;X4,X3]([*:1])[c]21")
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "Not a 1-benzopyran"

    # 2. Define an aryl group, allowing any aryl substituent
    aryl_pattern = Chem.MolFromSmarts("[cX3]1[cX3][cX3][cX3][cX3][cX3]1")

    # 3. Check for direct attachment at position 3
    # Create SMARTS pattern that explicitly connects the aryl to the C3 position
    # The C atom at position 3 is connected to any aryl via a single bond "-"
    combined_pattern = Chem.MolFromSmarts("c1ccccc2[O;X2][C;X4,X3](-[cX3]1[cX3][cX3][cX3][cX3][cX3]1)[c]21")
    
    if not mol.HasSubstructMatch(combined_pattern):
        return False, "No aryl group at position 3"

    return True, "1-Benzopyran with an aryl substituent at position 3"