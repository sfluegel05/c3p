"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies: lactol
Definition: Cyclic hemiacetals formed by intramolecular addition of a hydroxy group 
to an aldehydic or ketonic carbonyl group. In a lactol the carbonyl is converted into 
an sp3 “hemiacetal” carbon which typically carries one –OH and one –OR substituent 
(with the –OR being part of the ring).

This version relaxes the requirement on the candidate carbon’s hydrogen count (allowing 0 or 1)
and ensures that the two key oxygen neighbors are not the same.
"""

from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines whether the given SMILES represents a molecule featuring a lactol center.

    Heuristic approach:
      1. Parse the SMILES string and add explicit hydrogens.
      2. Loop over carbon atoms that are in a ring and sp3-hybridized.
      3. For each candidate carbon, require that:
           - It carries 0 or 1 hydrogen(s) (to cover both ketol and aldo cases).
           - It has at least 2 oxygen neighbors.
      4. Among its oxygen neighbors, we look for:
           - One oxygen (OH_candidate) that carries at least one explicit hydrogen (an –OH group).
           - A second oxygen (OR_candidate) that is in a ring (serving as the –OR of the cycle).
             (Ensure it’s a different atom than the OH_candidate.)
      5. Finally, verify that the candidate carbon and the OR_candidate appear together in
         at least one ring of size 5–7. 

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a lactol center is detected, False otherwise.
        str: A reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to reliably check for -OH groups and H count.
    mol = Chem.AddHs(mol)
    
    # Retrieve ring information: each ring as a tuple of atom indices.
    ring_info = mol.GetRingInfo().AtomRings()
    valid_ring_sizes = {5, 6, 7}  # typical ring sizes that can form lactols

    # Loop over all atoms to search for a candidate lactol center.
    for atom in mol.GetAtoms():
        # Consider only carbon atoms.
        if atom.GetAtomicNum() != 6:
            continue
        # Must be in a ring and sp3-hybridized.
        if not atom.IsInRing():
            continue
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        
        # For a hemiacetal, the carbon may have 0 or 1 hydrogen(s).
        num_H = atom.GetTotalNumHs()
        if num_H not in (0, 1):
            continue
        
        # Gather oxygen neighbors.
        oxy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
        if len(oxy_neighbors) < 2:
            continue
        
        # Look for two distinct oxygen roles:
        OH_candidate = None    # oxygen that carries an explicit hydrogen (-OH group)
        OR_candidate = None    # oxygen that is part of a ring (the -OR group)
        for oxy in oxy_neighbors:
            # Check if oxygen has at least one explicit hydrogen.
            oxy_has_H = any(nbr.GetAtomicNum() == 1 for nbr in oxy.GetNeighbors())
            # Check if oxygen is in a ring.
            oxy_in_ring = oxy.IsInRing()
            # Assign an -OH candidate if not already set.
            if oxy_has_H and OH_candidate is None:
                OH_candidate = oxy
                continue  # continue to check for OR_candidate
            # Assign an in-ring oxygen candidate if not already set.
            if oxy_in_ring and OR_candidate is None:
                OR_candidate = oxy
        # Ensure we have two distinct oxygens.
        if OH_candidate is None or OR_candidate is None or OH_candidate.GetIdx() == OR_candidate.GetIdx():
            continue

        # Now, verify that the candidate carbon and the OR_candidate share a ring of allowed size.
        c_idx = atom.GetIdx()
        o_idx = OR_candidate.GetIdx()
        shared_ring = False
        for ring in ring_info:
            if c_idx in ring and o_idx in ring and len(ring) in valid_ring_sizes:
                shared_ring = True
                break
        if not shared_ring:
            continue

        # If we have reached here, we have found a candidate lactol center.
        return True, f"Found cyclic hemiacetal (lactol) center at carbon atom index {c_idx}"
    
    return False, "No cyclic hemiacetal (lactol) center found"

# Optional: simple test when run as a script.
if __name__ == "__main__":
    # Test on beta-ascarylopyranose (a known lactol)
    test_smiles = "C[C@@H]1O[C@H](O)[C@H](O)C[C@H]1O"
    result, reason = is_lactol(test_smiles)
    print(result, reason)