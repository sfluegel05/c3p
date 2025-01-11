"""
Classifies: CHEBI:35358 sulfonamide
"""
from rdkit import Chem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide contains the functional group RS(=O)2NR'2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfonamide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the sulfonamide pattern
    sulfonamide_pattern = Chem.MolFromSmarts("S(=O)(=O)[N]")
    
    # Check for the sulfonamide group
    if not mol.HasSubstructMatch(sulfonamide_pattern):
        return False, "No sulfonamide group found"
    
    return True, "Contains a sulfonamide group"