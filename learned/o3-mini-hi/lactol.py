"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies: lactol 
Definition: Cyclic hemiacetals formed by intramolecular addition of a hydroxy group 
to an aldehydic or ketonic carbonyl group. In a lactol the carbonyl is converted into 
an sp3 “hemiacetal” carbon which typically carries one –OH and one –OR substituent (the –OR 
being part of the ring). This code inspects the molecule for a candidate lactol center.
"""

from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    
    Our heuristic approach is as follows:
      1. Parse the SMILES and add explicit hydrogens.
      2. Loop over every carbon atom that is (a) in a ring and (b) sp3-hybridized.
      3. For each such carbon, find oxygen neighbors. For a lactol center there
         should be exactly two oxygen neighbors:
           • One oxygen should be an –OH (having at least one hydrogen attached).
           • One oxygen should be an ether oxygen that is in the same small ring as the carbon.
      4. Check that the carbon and that “in‐ring” oxygen share a ring of size 5–7.
      5. If such a candidate is found, return True and the candidate’s atom index.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a lactol (cyclic hemiacetal) center is identified, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Make sure explicit hydrogens are added so we can detect –OH groups.
    mol = Chem.AddHs(mol)
    
    # Get ring information list (a tuple of atom index tuples for each ring)
    ring_info = mol.GetRingInfo().AtomRings()
    valid_ring_sizes = {5, 6, 7}  # typical small rings for lactols
    
    # Loop over all carbon atoms in the molecule.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # Consider only carbons in a ring and with sp3 hybridization.
        if not atom.IsInRing():
            continue
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
            
        # Gather oxygen neighbors of this carbon.
        oxy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
        # In a proper hemiacetal we expect exactly two oxygen neighbors.
        if len(oxy_neighbors) != 2:
            continue
        
        # Initialize flags and record candidate indices.
        oh_candidate = None   # oxygen that is –OH (has at least one hydrogen attached)
        or_candidate = None   # oxygen that is part of a ring with the carbon
        
        # Examine each oxygen neighbor.
        for oxy in oxy_neighbors:
            # Determine if the oxygen has any hydrogen neighbor (explicit H)
            has_h = any(nbr.GetAtomicNum() == 1 for nbr in oxy.GetNeighbors())
            # Check if this oxygen is in a ring.
            in_ring = oxy.IsInRing()
            # We assign based on the rule: one –OH (free hydroxyl, having an H) and one oxygen from the ring.
            if has_h and (oh_candidate is None):
                oh_candidate = oxy
            elif in_ring and (or_candidate is None):
                or_candidate = oxy
        
        # We need exactly one candidate for each role.
        if oh_candidate is None or or_candidate is None:
            continue
        
        # Confirm that the carbon and the ring oxygen actually belong to a shared ring
        # of size in the valid range.
        c_idx = atom.GetIdx()
        o_idx = or_candidate.GetIdx()
        in_valid_ring = False
        for ring in ring_info:
            if c_idx in ring and o_idx in ring and len(ring) in valid_ring_sizes:
                in_valid_ring = True
                break
        
        if not in_valid_ring:
            continue
        
        # At this point, we have found a candidate lactol center.
        return True, f"Found cyclic hemiacetal center at carbon atom index {c_idx}"
    
    return False, "No cyclic hemiacetal (lactol) center found"

# Example usage (optional):
if __name__ == "__main__":
    # Test on beta-ascarylopyranose (known lactol)
    smiles_example = "C[C@@H]1O[C@H](O)[C@H](O)C[C@H]1O"
    result, reason = is_lactol(smiles_example)
    print(result, reason)