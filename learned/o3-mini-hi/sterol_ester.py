"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: Sterol Ester
Definition: A steroid ester obtained by formal condensation of the carboxy group 
of any carboxylic acid with the 3-hydroxy group of a sterol.

Heuristic:
  1. The molecule must contain an ester functional group, i.e. a C(=O)O fragment.
  2. Among the ester groups, the oxygen that formerly belonged to the sterol’s
     3-hydroxy group (the ester oxygen not bound to the carbonyl) should be attached 
     to a ring atom.
  3. The molecule is required to have a steroid-type ring system. Many sterols contain 
     a cyclopentanoperhydrophenanthrene nucleus (four fused rings consisting of three 6-membered 
     rings and one 5-membered ring). We enforce that by requiring at least one 5-membered ring 
     and at least three 6-membered rings.
     
This heuristic does not cover all edge cases but should capture many examples of sterol ester structures.
"""
from rdkit import Chem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    
    A sterol ester is defined as a steroid ester formed by condensation of any 
    carboxylic acid with the 3-hydroxy group of a sterol. In other words, there 
    should be an ester (–C(=O)O–) where the oxygen that was formerly hydroxyl is attached
    to a steroid (fused four-ring) moiety.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is classified as a sterol ester, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, check for the presence of an ester group.
    # This SMARTS pattern looks for a C(=O)O fragment.
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester functional group (C(=O)O) found"
    
    # Next, check for a steroid nucleus.
    # A typical sterol has a cyclopentanoperhydrophenanthrene nucleus: four fused rings 
    # (usually three 6-membered rings and one 5-membered ring).
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    num_5_membered = sum(1 for ring in rings if len(ring) == 5)
    num_6_membered = sum(1 for ring in rings if len(ring) == 6)
    if num_5_membered < 1 or num_6_membered < 3:
        return False, ("Steroid nucleus not detected – "
                       "expected at least one 5-membered ring and three 6-membered rings, "
                       f"found {num_5_membered} and {num_6_membered}, respectively")
    
    # Among the ester groups, look for one that is likely formed by the condensation 
    # at the 3-hydroxy of a sterol.
    # In an ester, the oxygen attached to the carbonyl (C(=O)O) will have two neighbors.
    # One neighbor (the carbonyl C) belongs to the acid part; the other,
    # which should be the former sterol hydroxyl group, should be attached to a ring.
    sterol_ester_found = False
    for match in ester_matches:
        # The SMARTS "C(=O)O" should match at least two atoms. 
        # We assume here:
        #   match[0] is the carbonyl carbon
        #   match[1] is the oxygen (the ester oxygen, coming from the alcohol)
        if len(match) < 2:
            continue
        carbonyl_idx = match[0]
        ester_oxygen_idx = match[1]
        ester_oxygen = mol.GetAtomWithIdx(ester_oxygen_idx)
        
        # Get neighbors of this oxygen atom. It should be connected to the carbonyl carbon
        # and one additional atom. The additional atom should be part of a ring,
        # presumably coming from the sterol moiety.
        neighbors = [nbr.GetIdx() for nbr in ester_oxygen.GetNeighbors()]
        if len(neighbors) != 2:
            continue  # this ester group is unusual; skip it
        
        # Identify the neighbor that is not the carbonyl carbon.
        if neighbors[0] == carbonyl_idx:
            steroid_attachment = neighbors[1]
        elif neighbors[1] == carbonyl_idx:
            steroid_attachment = neighbors[0]
        else:
            continue
        
        # Check if the attachment point is part of a ring.
        attachment_atom = mol.GetAtomWithIdx(steroid_attachment)
        if attachment_atom.IsInRing():
            sterol_ester_found = True
            break
    
    if not sterol_ester_found:
        return False, ("Ester group found but the oxygen from the ester does not appear to "
                       "connect to a ring system (expected for a sterol-derived alcohol)")
    
    return True, "Contains an ester group and a ring system consistent with a sterol nucleus"

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided example SMILES
    smiles_example = "O([C@@H]1CC=2[C@@](C3C(C4[C@@](C(CC4)[C@@H](CCCC(C)C)C)(CC3)C)CC2)(CC1)C)C(CCCCCCCCC/C=C\\C/C=C\\CCCCC)=O"
    result, reason = is_sterol_ester(smiles_example)
    print("Result:", result)
    print("Reason:", reason)