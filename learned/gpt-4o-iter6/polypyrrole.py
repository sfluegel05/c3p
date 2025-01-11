"""
Classifies: CHEBI:38077 polypyrrole
"""
from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is defined as a compound composed of two or more pyrrole units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypyrrole, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a pyrrole unit
    pyrrole_pattern = Chem.MolFromSmarts("n1cccc1")
    if pyrrole_pattern is None:
        return False, "Failed to generate SMARTS pattern for pyrrole"

    # Find all pyrrole matches in the molecule
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    
    # Check if there are two or more pyrrole units
    num_pyrroles = len(pyrrole_matches)
    if num_pyrroles >= 2:
        return True, f"Contains {num_pyrroles} pyrrole units"
    else:
        return False, "Contains less than two pyrrole units"

# Example usage:
smiles_example = "C1=CNC=C1C1=CC=CN=C1"
is_polypyrrole(smiles_example)