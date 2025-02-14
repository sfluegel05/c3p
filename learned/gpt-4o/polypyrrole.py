"""
Classifies: CHEBI:38077 polypyrrole
"""
from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is characterized by two or more pyrrole units,
    which are five-membered rings with nitrogen.

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
    
    # Define the pyrrole pattern
    pyrrole_pattern = Chem.MolFromSmarts('[nH]1cccc1')
    if not pyrrole_pattern:
        return None, "Invalid pyrrole SMARTS pattern"
    
    # Find matches for pyrrole units
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    
    # Check if there are two or more distinct pyrrole units.
    # Avoid multiple matches for the same pyrrole by using unique atom indices.
    unique_matches = set(tuple(sorted(match)) for match in pyrrole_matches)
    
    if len(unique_matches) >= 2:
        return True, f"Contains {len(unique_matches)} pyrrole units"

    return False, "Less than 2 pyrrole units found"