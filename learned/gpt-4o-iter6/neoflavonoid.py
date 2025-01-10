"""
Classifies: CHEBI:71971 neoflavonoid
"""
from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is any 1-benzopyran with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for 1-benzopyran core with aryl group at position 4
    benzopyran_aryl_pattern = Chem.MolFromSmarts("c1cc2occcc2c1-c3ccccc3")  # 1-benzopyran with an aryl group at position 4

    # Check for the presence of the specific 1-benzopyran with aryl pattern
    if mol.HasSubstructMatch(benzopyran_aryl_pattern):
        return True, "Contains 1-benzopyran with an aryl substituent at position 4"

    return False, "Does not match 1-benzopyran with an aryl substituent at position 4"