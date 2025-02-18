"""
Classifies: CHEBI:38077 polypyrrole
"""
"""
Classifies: CHEBI:xxxxx polypyrrole
"""
from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole contains two or more pyrrole units (5-membered aromatic rings with one nitrogen).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypyrrole, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define SMARTS pattern for a pyrrole ring: 5-membered aromatic ring with exactly one nitrogen
    pyrrole_smarts = '[n;ar]1[c;ar][c;ar][c;ar][c;ar]1'
    pyrrole_pattern = Chem.MolFromSmarts(pyrrole_smarts)
    
    if not pyrrole_pattern:
        return False, "Failed to initialize pyrrole pattern"
    
    # Find all unique pyrrole rings (avoid overlapping matches)
    matches = mol.GetSubstructMatches(pyrrole_pattern)
    unique_rings = set()
    for match in matches:
        # Sort atom indices to handle different atom orders in SMARTS matches
        unique_rings.add(tuple(sorted(match)))
    
    pyrrole_count = len(unique_rings)
    
    if pyrrole_count >= 2:
        return True, f"Contains {pyrrole_count} pyrrole rings"
    else:
        return False, f"Found {pyrrole_count} pyrrole rings, need at least 2"