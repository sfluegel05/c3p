"""
Classifies: CHEBI:38077 polypyrrole
"""
from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is a compound composed of two or more pyrrole units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypyrrole, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None  # Invalid SMILES, return None, None
    
    # Define the pyrrole substructure using SMARTS
    pyrrole_smarts = "[nX1]1[cX4][cX4][cX4][cX4]1" # Modified SMARTS to find a pyrrole ring
    pyrrole_pattern = Chem.MolFromSmarts(pyrrole_smarts)

    # Count the number of pyrrole rings
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    num_pyrroles = len(pyrrole_matches)

    # Classify based on the count of pyrrole rings
    if num_pyrroles >= 2:
         return True, f"Contains {num_pyrroles} pyrrole units, therefore a polypyrrole."
    else:
        return False, f"Contains only {num_pyrroles} pyrrole unit(s), therefore not a polypyrrole."