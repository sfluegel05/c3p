"""
Classifies: CHEBI:30527 flavin
"""
from rdkit import Chem

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin based on its SMILES string.
    A flavin is a derivative of dimethylisoalloxazine with a substituent on the 10th position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavin, False otherwise
        str: Reason for classification
    """
    # Convert SMILES string to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for dimethylisoalloxazine core
    # This pattern should ideally identify the heterocyclic core correctly.
    core_pattern = Chem.MolFromSmarts("Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C)c2cc1C")
    
    # Check for core structure
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core flavin structure (dimethylisoalloxazine) missing"
    
    # Check for substituent at 10th position (attached to the N atom in the core)
    substituent_pattern = Chem.MolFromSmarts("Cc1cc2nc3[nX3]c(=O)n(c3=O)c2cc1C")
    if not mol.HasSubstructMatch(substituent_pattern):
        return False, "No substituent detected at the 10th position"
    
    return True, "Valid flavin structure detected"