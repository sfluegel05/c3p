"""
Classifies: CHEBI:38077 polypyrrole
"""
from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is defined as a compound with two or more pyrrole units.

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

    # Define pyrrole pattern (5-membered aromatic ring with NH)
    pyrrole_pattern = Chem.MolFromSmarts("[n&H1]1cccc1")  # Ensure nitrogen is part of an aromatic system
    
    # Try multiple pyrrole-like patterns if necessary
    pyrrole_alt_pattern = Chem.MolFromSmarts("[n]1ccc[cH0]1")  # For fused systems or aromatic checks

    # Find matches for pyrrole units
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    alt_matches = mol.GetSubstructMatches(pyrrole_alt_pattern)
    
    # Calculate the total unique matches
    unique_matches = set(pyrrole_matches) | set(alt_matches)

    # Check if there are two or more pyrrole units
    if len(unique_matches) >= 2:
        return True, f"Contains {len(unique_matches)} pyrrole units"

    return False, "Less than 2 pyrrole units found"