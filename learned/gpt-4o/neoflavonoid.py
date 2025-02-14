"""
Classifies: CHEBI:71971 neoflavonoid
"""
from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is defined as any 1-benzopyran with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1-benzopyran core pattern (chromene structure)
    benzopyran_pattern = Chem.MolFromSmarts('c1ccc2occ(C)c2c1')  # basic chromene structure
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "1-benzopyran core not found"

    # Aryl substituent pattern at position 4
    aryl_pattern = Chem.MolFromSmarts('c1ccc(cc1)-c2ccc3occc3c2')  # flexible aryl pattern
    if not mol.HasSubstructMatch(aryl_pattern):
        return False, "Aryl group not found at position 4"

    return True, "Contains a 1-benzopyran core with an aryl group at position 4"