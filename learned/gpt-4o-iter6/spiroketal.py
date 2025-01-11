"""
Classifies: CHEBI:72600 spiroketal
"""
from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal is defined as a spiro center where a carbon is the only shared atom between two rings,
    each containing at least one oxygen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a spiroketal, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for spiro center
    # Look for a carbon spiro center that is connected to two distinct rings
    # Each ring should contain an oxygen atom
    spiroketal_smarts = "[C]([O])([*])"
    spiroketal_pattern = Chem.MolFromSmarts(spiroketal_smarts)
    if not spiroketal_pattern:
        return None, "Invalid SMARTS pattern"

    # Check for pattern in the molecule
    matches = mol.GetSubstructMatches(spiroketal_pattern)
    ring_info = mol.GetRingInfo()

    for match in matches:
        carbon_idx, o1_idx, _ = match
        
        # Ensure the carbon atom is part of exactly two rings
        rings_with_carbon = [ring for ring in ring_info.AtomRings() if carbon_idx in ring]
        if len(rings_with_carbon) != 2:
            continue
        
        # Verify that the oxygen atom is part of only one of the two rings
        o1_rings = [ring for ring in rings_with_carbon if o1_idx in ring]
        if len(o1_rings) != 1:
            continue
        
        return True, "Found a spiroketal with a carbon connecting two distinct ring oxygens"

    return False, "No spiroketal structure found"