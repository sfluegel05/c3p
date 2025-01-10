"""
Classifies: CHEBI:38077 polypyrrole
"""
"""
Classifies: polypyrrole
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

    # Define pyrrole pattern
    pyrrole_pattern = Chem.MolFromSmarts("n1cccc1")
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    
    # Count number of pyrrole substructures
    num_pyrrole_units = len(pyrrole_matches)
    
    if num_pyrrole_units < 2:
        return False, f"Contains only {num_pyrrole_units} pyrrole units, need 2 or more for polypyrrole"

    return True, f"Contains {num_pyrrole_units} pyrrole units, qualifies as polypyrrole"