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
    
    # Refining the SMARTS pattern for spiroketal
    # SMARTS for spiroketal requires a carbon spiro center with two oxygens, each in a different ring
    spiroketal_smarts = "[C]([O][R])[O][R]"  # Improved SMARTS, see if a single carbon joins two ether oxygens of different rings
    spiroketal_pattern = Chem.MolFromSmarts(spiroketal_smarts)
    if not spiroketal_pattern:
        return None, "Invalid SMARTS pattern"
    
    # Find matches
    matches = mol.GetSubstructMatches(spiroketal_pattern)
    ring_info = mol.GetRingInfo()
    
    for match in matches:
        # The first match index is the spiro carbon, followed by two oxygens
        carbon_idx = match[0]
        oxygen_indices = match[1:3]
        
        # Check that each oxygen is part of a distinct ring structure
        o_ring_assignments = [ring_info.AtomRingsFromIdx(o_idx) for o_idx in oxygen_indices]
        
        # Ensure oxygens are in different rings
        if all(len(oxygen_rings) > 0 for oxygen_rings in o_ring_assignments):
            first_ring_set, second_ring_set = [set(ring) for rings in o_ring_assignments for ring in rings]
            
            # Confirm they are in separate ring systems
            if len(set.intersection(first_ring_set, second_ring_set)) == 0:
                return True, "Found a spiroketal with distinct ring linkages via a spiro center"
            
    return False, "No spiroketal structure found"