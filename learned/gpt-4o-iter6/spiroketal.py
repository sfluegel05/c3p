"""
Classifies: CHEBI:72600 spiroketal
"""
from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal contains a spiro center that is the only common atom of two rings,
    with each ring containing one of the connected oxygen atoms.

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
    
    # Define SMARTS pattern for spiroketal
    # This SMARTS looks for a carbon spiro center with two oxygen atoms each in different rings
    spiroketal_smarts = "[C]([O])([O])"
    spiroketal_pattern = Chem.MolFromSmarts(spiroketal_smarts)
    if not spiroketal_pattern:
        return None, "Invalid SMARTS pattern"

    # Find matches
    matches = mol.GetSubstructMatches(spiroketal_pattern)
    ring_info = mol.GetRingInfo()
    
    for match in matches:
        # The first match index is the spiro carbon, followed by two oxygens
        carbon_idx, o1_idx, o2_idx = match
        
        # Get rings each oxygen is part of
        oxygen_rings_1 = set(ring for ring in ring_info.AtomRings() if o1_idx in ring)
        oxygen_rings_2 = set(ring for ring in ring_info.AtomRings() if o2_idx in ring)
        
        # Check that each oxygen is in a different ring
        if oxygen_rings_1 and oxygen_rings_2 and oxygen_rings_1 != oxygen_rings_2:
            return True, "Found a spiroketal with a carbon connecting two distinct ring oxygens"
    
    return False, "No spiroketal structure found"