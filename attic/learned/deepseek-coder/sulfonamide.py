"""
Classifies: CHEBI:35358 sulfonamide
"""
"""
Classifies: CHEBI:35358 sulfonamide
"""
from rdkit import Chem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide is characterized by the functional group RS(=O)2NR'2.

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

    # Define the sulfonamide pattern: S(=O)(=O)-N
    sulfonamide_pattern = Chem.MolFromSmarts("[SX4](=[OX1])(=[OX1])[NX3]")
    
    # Check if the molecule contains the sulfonamide pattern
    if mol.HasSubstructMatch(sulfonamide_pattern):
        return True, "Contains the sulfonamide functional group RS(=O)2NR'2"
    else:
        return False, "No sulfonamide functional group found"