"""
Classifies: CHEBI:72600 spiroketal
"""
#!/usr/bin/env python
"""
Classifies: spiroketal
Definition: A cyclic ketal in which the ketal carbon is the only common atom of two rings.
A true spiroketal features a carbon (ideally sp³ but we now relax that requirement)
that is tetra‐bonded to exactly two oxygen and two carbon heavy atoms. Moreover,
there must be two rings (of a reasonable minimum size) that each contain the candidate carbon
and exactly one of its oxygen neighbors so that the rings share only the candidate atom.
"""

from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    
    The algorithm uses the following criteria:
      1. The candidate spiro center must be a carbon (atomic number 6).
      2. It must have exactly four heavy-atom neighbors (ignoring implicit hydrogens),
         exactly two of which are oxygens (atomic number 8) and two carbons.
      3. The candidate must belong to at least two rings (ignoring rings with less than 4 atoms).
      4. For each ring that contains the candidate, we count how many of the candidate’s oxygen
         neighbors are present. For a true spiroketal each of the two rings that come together at
         the candidate should contain exactly one (and a different) oxygen neighbor.
      5. Finally, we require that the pair of rings chosen share only the candidate atom.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a spiroketal, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information as a list of sets of atom indices.
    rings = [set(r) for r in mol.GetRingInfo().AtomRings() if len(r) >= 4]
    if not rings:
        return False, "No rings detected in the molecule (or no rings of size>=4)"
    
    # Loop over atoms searching for candidate spiro centers.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # (Relax hybridization requirement so that stereochemically annotated atoms are not skipped.)
    
        # Use the explicit neighbors (heavy atoms only)
        neighbors = atom.GetNeighbors()
        if len(neighbors) != 4:
            continue  # Must be tetravalent when considering heavy atoms.
    
        # Count oxygen and carbon neighbors.
        oxygen_neighbors = []
        carbon_neighbors = []
        for nbr in neighbors:
            anum = nbr.GetAtomicNum()
            if anum == 8:
                oxygen_neighbors.append(nbr.GetIdx())
            elif anum == 6:
                carbon_neighbors.append(nbr.GetIdx())
        if len(oxygen_neighbors) != 2 or len(carbon_neighbors) != 2:
            continue  # Wrong substitution pattern.
    
        cand_idx = atom.GetIdx()
        # Gather rings that contain the candidate.
        candidate_rings = [r for r in rings if cand_idx in r]
        if len(candidate_rings) < 2:
            continue
    
        # For each candidate ring, count how many oxygen neighbors occur in that ring.
        # We want rings that (a) contain exactly one of the candidate’s oxygen neighbors.
        ring_to_ox = {}  # Map ring (as frozenset) -> oxygen neighbor index present in that ring.
        for r in candidate_rings:
            # Find intersection between the ring and the oxygen neighbor set.
            ox_in_ring = set(oxygen_neighbors).intersection(r)
            if len(ox_in_ring) == 1:
                # Save the ring (using frozenset as key) along with the oxygen index present.
                ring_to_ox[frozenset(r)] = list(ox_in_ring)[0]
        
        # To be a spiroketal we need two distinct rings that have different oxygen neighbors
        found_pair = None
        ring_list = list(ring_to_ox.items())
        for i in range(len(ring_list)):
            for j in range(i+1, len(ring_list)):
                ring1, ox1 = ring_list[i]
                ring2, ox2 = ring_list[j]
                if ox1 == ox2:
                    # Both rings have the same oxygen from the candidate, so they do not
                    # split the two distinct oxygen substituents.
                    continue
                # Check that the two rings share only the candidate atom.
                if ring1.intersection(ring2) == {cand_idx}:
                    found_pair = (ring1, ring2, ox1, ox2)
                    break
            if found_pair:
                break
        if found_pair:
            ring1, ring2, ox1, ox2 = found_pair
            # Identify the carbon neighbors for reporting (order arbitrary)
            carbons = carbon_neighbors
            return True, (f"Found spiroketal center at atom index {cand_idx} with oxygen neighbors "
                          f"{oxygen_neighbors[0]} and {oxygen_neighbors[1]} and carbon neighbors "
                          f"{carbons[0]} and {carbons[1]}. Ring one (size {len(ring1)}) contains oxygen {ox1}; "
                          f"ring two (size {len(ring2)}) contains oxygen {ox2}; rings share only candidate atom.")
    
    return False, "No spiroketal pattern found in the molecule"

# Example usage:
if __name__ == "__main__":
    # Test with one known positive (chrysopine) and a known false positive example.
    test_smiles_pos = "NC(=O)CC[C@@H]1NCC2(OC[C@@H](O)[C@@H](O)[C@@H]2O)OC1=O"
    result, reason = is_spiroketal(test_smiles_pos)
    print("Positive test:", result, reason)
    
    test_smiles_neg = "O=C1O[C@@]2(OC[C@]3([C@@H]2O)[C@@H]4C(=C)CC[C@H]3C4)C=C1C"  # Example false positive
    result2, reason2 = is_spiroketal(test_smiles_neg)
    print("Negative test:", result2, reason2)