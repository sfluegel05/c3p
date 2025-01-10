"""
Classifies: CHEBI:72600 spiroketal
"""
from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal contains a spiro center that is the only common atom of two rings,
    with two attached oxygen atoms each part of a different ring.

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
    
    # Define a SMARTS pattern for spiroketal
    # The pattern can be defined as a carbon atom involved as a spiro center having
    # two connections to oxygens with those oxygens being in separate rings
    spiroketal_smarts = "[C]1([O])([O])" # Improvised SMARTS (this likely needs refinement)
    
    spiroketal_pattern = Chem.MolFromSmarts(spiroketal_smarts)
    if not spiroketal_pattern:
        return None, "Invalid SMARTS pattern"
    
    # Find matches
    matches = mol.GetSubstructMatches(spiroketal_pattern)
    # Ring information
    ring_info = mol.GetRingInfo()
    
    for match in matches:
        # Check ring count for oxygen attachments
        # Assuming match provides indices, needs to verify each O in different ring structs
        carbon_idx = match[0]
        o_indices = match[1:3]
        if all(ring_info.NumAtomRings(o_idx) >= 1 for o_idx in o_indices):
            o_rings = [ring_info.AtomRingsFromIdx(o_idx) for o_idx in o_indices]
            if all(len(o_ring_set) > 0 for o_ring_set in o_rings):
                # Ensure distinct ring sets
                ring_sets = [set(ring) for rings in o_rings for ring in rings]
                # Check that each oxygen is part of a separate ring
                if len(set.intersection(*ring_sets)) == 0:
                    return True, "Found a spiroketal with the correct linkage"
            
    return False, "No spiroketal found"