"""
Classifies: CHEBI:35358 sulfonamide
"""
from rdkit import Chem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfonamide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the sulfonamide SMARTS pattern
    sulfonamide_pattern = Chem.MolFromSmarts("S(=O)(=O)[N]")

    # Check if the molecule has the sulfonamide substructure
    if mol.HasSubstructMatch(sulfonamide_pattern):
        return True, "Contains sulfonamide group"

    return False, "No sulfonamide group found"