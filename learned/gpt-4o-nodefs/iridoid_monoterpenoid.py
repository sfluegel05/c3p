"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an iridoid monoterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a cyclopentanopyran system
    # Note: this is simplified and may not cover all varieties of iridoid monoterpenoids
    cyclopentanopyran_smarts = "C1CCC2C(C1)OC=C2"  # Example of pattern for a cyclopentanopyran
    
    # Convert the SMARTS to a molecule object
    cyclopentanopyran_pattern = Chem.MolFromSmarts(cyclopentanopyran_smarts)
    if mol.HasSubstructMatch(cyclopentanopyran_pattern):
        return True, "Compound contains a cyclopentanopyran framework indicative of iridoid monoterpenoids"
    
    return False, "No cyclopentanopyran framework found; does not match typical iridoid monoterpenoid structure"