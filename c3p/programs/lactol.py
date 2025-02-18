"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies: lactol
Definition: Cyclic hemiacetals formed by intramolecular addition of a hydroxy group 
to an aldehydic or ketonic carbonyl group. In a lactol the carbonyl is converted into 
an sp3 “hemiacetal” carbon which typically carries one –OH and one –OR substituent 
(with the –OR being part of the ring). 

This improved version uses additional checks:
  • Adds explicit hydrogens.
  • Loops over carbons that are in a ring, sp3‐hybridized, and that carry exactly 1 hydrogen 
    (as expected for a hemiacetal carbon, CH(OH)(OR)).
  • Requires that the candidate carbon has at least one oxygen neighbor that has an explicit 
    hydrogen (an –OH) and at least one oxygen neighbor that is in a ring (the –OR of the cyclic structure).
  • Confirms that the candidate carbon and that “in‐ring” oxygen occur together in a ring 
    of size 5–7.
If such a candidate is found, returns True and the candidate’s atom index.
If no candidate is found (or the SMILES can’t be parsed) returns False.
"""

from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if the given SMILES string represents a molecule featuring a lactol center.
    
    Our heuristic approach:
      1. Parse the SMILES string and add explicit hydrogens.
      2. For each carbon atom:
           - It must be in a ring, sp3-hybridized,
           - It should have exactly one hydrogen attached (as expected for a hemiacetal CH),
             and a total of at least 2 non-hydrogen neighbors.
      3. Among its neighbors, we require:
           - At least one oxygen neighbor with an explicit hydrogen (–OH group).
           - At least one oxygen neighbor that is itself in a ring (acting as the ring –OR).
      4. Verify that the candidate carbon and the ring oxygen share at least one ring 
         of size 5–7.
      5. If a candidate is found, return True along with a message showing the candidate’s index.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a lactol (cyclic hemiacetal) center is detected, False otherwise.
        str: A reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens needed to check for -OH groups and hydrogen count on carbon.
    mol = Chem.AddHs(mol)
    
    # Retrieve ring information: each ring as a tuple of atom indices.
    ring_info = mol.GetRingInfo().AtomRings()
    valid_ring_sizes = {5, 6, 7}  # typical ring sizes that can form lactols
    
    # Loop over all atoms to select candidate lactol centers.
    for atom in mol.GetAtoms():
        # Consider only carbon atoms.
        if atom.GetAtomicNum() != 6:
            continue
        # Must be in a ring and sp3-hybridized.
        if not atom.IsInRing():
            continue
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        # Check that the carbon has exactly 1 hydrogen attached.
        # (For a hemiacetal formed from a carbonyl, we expect CH(OH)(OR) pattern.)
        num_H = atom.GetTotalNumHs()
        if num_H != 1:
            continue
        
        # Look at the atom's neighbors; we want to find oxygen neighbors.
        oxy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
        # (We do not insist on exactly two oxygens because sometimes extra substituents exist,
        # but we require at least two and that the roles are filled.)
        if len(oxy_neighbors) < 2:
            continue
        
        # Initialize indicators for the two roles.
        oh_found = False       # oxygen with at least one explicit hydrogen (-OH)
        ring_oxygen_found = None   # oxygen that is in a ring (likely the -OR group)
        
        for oxy in oxy_neighbors:
            # Determine if oxygen bears at least one explicit hydrogen.
            oxy_has_H = any(nbr.GetAtomicNum() == 1 for nbr in oxy.GetNeighbors())
            # Determine if oxygen is in a ring.
            oxy_in_ring = oxy.IsInRing()
            # We expect one oxygen to be an –OH:
            if oxy_has_H and not oh_found:
                oh_found = True
                continue  # check next neighbor
            # And we expect at least one oxygen to be part of a ring.
            if oxy_in_ring and ring_oxygen_found is None:
                ring_oxygen_found = oxy
        
        # If we did not meet both roles, skip this carbon.
        if not oh_found or ring_oxygen_found is None:
            continue
        
        # Confirm that this carbon and the ring oxygen are in a shared ring of valid size.
        c_idx = atom.GetIdx()
        o_idx = ring_oxygen_found.GetIdx()
        in_valid_ring = False
        for ring in ring_info:
            if c_idx in ring and o_idx in ring and len(ring) in valid_ring_sizes:
                in_valid_ring = True
                break
        if not in_valid_ring:
            continue
        
        # Found a candidate lactol center.
        return True, f"Found cyclic hemiacetal (lactol) center at carbon atom index {c_idx}"
    
    return False, "No cyclic hemiacetal (lactol) center found"

# Optional: simple test when run as a script.
if __name__ == "__main__":
    # Test on beta-ascarylopyranose (a known lactol)
    test_smiles = "C[C@@H]1O[C@H](O)[C@H](O)C[C@H]1O"
    result, reason = is_lactol(test_smiles)
    print(result, reason)