"""
Classifies: CHEBI:72600 spiroketal
"""
#!/usr/bin/env python
"""
Classifies: spiroketal
Definition: A cyclic ketal in which the ketal carbon is the only common atom of two rings.
A true spiroketal features a tetra‐substituted carbon (the candidate spiro center)
that is bonded by single bonds to exactly two oxygens and two carbons, and that lies at the
junction of two rings (each of reasonable size) that each include one – and a different – oxygen neighbor.
"""

from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    
    The algorithm uses the following criteria:
      1. The candidate spiro center must be a carbon (atomic number 6) that is fully substituted (degree including implicit Hs equals 4)
         and must not have any implicit hydrogens.
      2. It must have exactly four explicit (heavy-atom) bonds: two to oxygen atoms (via SINGLE bonds) and two to carbon atoms.
      3. The candidate must be present in at least two rings (only rings of size >=5 are considered).
      4. For each ring containing the candidate, we count how many of the candidate’s oxygen neighbors occur in that ring.
         A valid ring should feature exactly one oxygen neighbor.
      5. Two rings that share only the candidate *and* that each include a different oxygen neighbor establish a spiroketal motif.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a spiroketal, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information as a list of sets of atom indices (only rings with at least 5 atoms).
    rings = [set(r) for r in mol.GetRingInfo().AtomRings() if len(r) >= 5]
    if not rings:
        return False, "No rings detected in the molecule (or no rings of size>=5)"

    # Loop over atoms looking for candidate spiro centers.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue

        # For a true ketal center, total degree (explicit + implicit) should be 4,
        # and there should be no implicit hydrogens attached.
        if atom.GetTotalDegree() != 4 or atom.GetNumImplicitHs() != 0:
            continue

        # Get explicit (heavy atom) neighbors.
        neighbors = atom.GetNeighbors()
        # We require exactly 4 explicit neighbors.
        if len(neighbors) != 4:
            continue

        # Check bonds to each neighbor: for oxygen neighbors, ensure the bond is single.
        oxygen_neighbors = []
        carbon_neighbors = []
        valid_candidate = True
        for nbr in neighbors:
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            anum = nbr.GetAtomicNum()
            if anum == 8:
                # For a valid ketal, the oxygen must be linked by a SINGLE bond.
                if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                    valid_candidate = False
                    break
                oxygen_neighbors.append(nbr.GetIdx())
            elif anum == 6:
                carbon_neighbors.append(nbr.GetIdx())
            else:
                # If any other heavy atom is present, skip.
                valid_candidate = False
                break
        if not valid_candidate:
            continue
        if len(oxygen_neighbors) != 2 or len(carbon_neighbors) != 2:
            continue  # Wrong substitution pattern for a ketal center.

        cand_idx = atom.GetIdx()
        # Gather rings that contain the candidate.
        candidate_rings = [r for r in rings if cand_idx in r]
        if len(candidate_rings) < 2:
            continue

        # For each candidate ring, count how many oxygen-neighbors (attached by a single bond) occur in that ring.
        # We want rings that contain exactly one of the candidate’s oxygen neighbors.
        ring_to_ox = {}  # Map: ring (as frozenset) -> oxygen neighbor index present in that ring.
        for r in candidate_rings:
            ox_in_ring = set(oxygen_neighbors).intersection(r)
            if len(ox_in_ring) == 1:
                # Save the ring (as frozenset) along with the oxygen index present.
                ring_to_ox[frozenset(r)] = list(ox_in_ring)[0]

        # Look for a pair of rings that use different oxygen neighbors and share only the candidate.
        ring_list = list(ring_to_ox.items())
        for i in range(len(ring_list)):
            for j in range(i + 1, len(ring_list)):
                ring1, ox1 = ring_list[i]
                ring2, ox2 = ring_list[j]
                if ox1 == ox2:
                    # Both rings use the same oxygen; not acceptable.
                    continue
                # Verify that the two rings share only the candidate atom.
                if ring1.intersection(ring2) == {cand_idx}:
                    # Found a valid pair.
                    carbons = carbon_neighbors
                    return True, (f"Found spiroketal center at atom index {cand_idx} with oxygen neighbors "
                                  f"{oxygen_neighbors[0]} and {oxygen_neighbors[1]} and carbon neighbors "
                                  f"{carbons[0]} and {carbons[1]}. Ring one (size {len(ring1)}) contains oxygen {ox1}; "
                                  f"ring two (size {len(ring2)}) contains oxygen {ox2}; rings share only candidate atom.")

    return False, "No spiroketal pattern found in the molecule"

# Example usage:
if __name__ == "__main__":
    # Positive test: chrysopine (should be classified as spiroketal).
    test_smiles_pos = "NC(=O)CC[C@@H]1NCC2(OC[C@@H](O)[C@@H](O)[C@@H]2O)OC1=O"
    res, reason = is_spiroketal(test_smiles_pos)
    print("Positive test:", res, reason)
    
    # A false positive example that was previously misclassified.
    test_smiles_neg = "O=C1O[C@@]2(OC[C@]3([C@@H]2O)[C@@H]4C(=C)CC[C@H]3C4)C=C1C"
    res2, reason2 = is_spiroketal(test_smiles_neg)
    print("Negative test:", res2, reason2)