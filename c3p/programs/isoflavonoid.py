"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: Any 1-benzopyran with an aryl substituent at position 3 (isoflavonoids)

An isoflavonoid is defined as a benzopyran (chromene) system fused with a benzene ring,
with one of the positions on the pyran (here, our candidate “position 3”) bearing an exocyclic aryl substituent.
Since isoflavonoids occur with a wide structural diversity, instead of a single SMARTS query
we use an algorithmic approach. We:
  1. Collect all 6-membered aromatic rings that contain exactly one oxygen (candidate “pyran” rings).
  2. For each candidate, check if it is fused (shares at least 2 atoms) with a 6-membered aromatic ring
     composed entirely of carbons (a benzene).
  3. Then, among the atoms of the candidate pyran that are not in the fused region,
     confirm that at least one carries a substituent that is part of a separate 6-membered aromatic ring.
If all these checks pass we return True.
"""

from rdkit import Chem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    
    An isoflavonoid is defined as any 1-benzopyran with an aryl substituent at position 3.
    This implementation first extracts rings from the molecule and then applies several filters:
      1. Find a candidate 6-membered aromatic ring that has exactly one oxygen.
      2. Check that the candidate is fused (at least 2 atoms) with a 6-membered aromatic ring that is
         composed solely of carbons (the benzene ring).
      3. For one or more atoms of the candidate (in particular, a position corresponding to C-3) that
         are not shared with the fused benzene, require that an external substituent is part of an aromatic ring.
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an isoflavonoid, False otherwise.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information (each ring is a tuple of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found"

    # Helper: filter rings that are 6-membered and fully aromatic.
    candidate_rings = []
    for ring in ring_info:
        if len(ring) == 6:
            # Check that every atom in the ring is aromatic
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                candidate_rings.append(set(ring))
    
    if not candidate_rings:
        return False, "No 6-membered aromatic rings found"
    
    # Function: count number of oxygen atoms within a set of atom indices
    def count_oxygen(ring_set):
        count = 0
        for idx in ring_set:
            if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8:
                count += 1
        return count

    # Now look for a candidate pyran ring: a 6-membered aromatic ring with exactly one oxygen.
    pyran_candidates = [ring for ring in candidate_rings if count_oxygen(ring) == 1]
    if not pyran_candidates:
        return False, "No 6-membered aromatic ring with exactly one oxygen (pyran candidate) found"
    
    # For each pyran candidate, look for a fused benzene ring (a 6-membered fully aromatic ring with only carbons)
    for pyran in pyran_candidates:
        for other in candidate_rings:
            # Skip if they are identical
            if pyran == other:
                continue
            # Check if other ring is purely carbons.
            if any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 for idx in other):
                continue
            # Determine the fusion: the intersection should have at least 2 atoms.
            fusion = pyran.intersection(other)
            if len(fusion) < 2:
                continue
            
            # Now we have a fused pair: pyran (chromene part) and benzene.
            fused_indices = pyran.union(other)
            # Now, in the pyran ring, consider atoms that are NOT in the fused intersection.
            nonfused_atoms = pyran.difference(fusion)
            # We expect that one of these atoms (often position 3 in a chromene) carries an external substituent 
            # that is an aryl group (i.e. belongs to a separate aromatic 6-membered ring).
            for idx in nonfused_atoms:
                atom = mol.GetAtomWithIdx(idx)
                # Look at neighbors (outside the pyran candidate)
                for neighbor in atom.GetNeighbors():
                    n_idx = neighbor.GetIdx()
                    if n_idx in pyran:
                        continue  # skip atoms inside the pyran
                    # Check if neighbor is aromatic carbon
                    if not (neighbor.GetIsAromatic() and neighbor.GetAtomicNum() == 6):
                        continue
                    # Now check if this neighbor is in some 6-membered aromatic ring that is disjoint from the entire fused system.
                    neighbor_in_aryl = False
                    for ring in candidate_rings:
                        # We want an aromatic ring that is not a subset of our fused system.
                        # In practice the aryl substituent is typically a benzene ring that has little overlap.
                        if n_idx in ring and len(ring) == 6 and ring.isdisjoint(fused_indices):
                            neighbor_in_aryl = True
                            break
                    if neighbor_in_aryl:
                        # We found an external aryl substituent on one of the non-fused positions.
                        return True, ("Molecule contains a benzopyran (chromene) core fused with a benzene ring and "
                                      "an external aryl substituent on a non-fused position, consistent with an isoflavonoid scaffold.")
    return False, "Scaffold not recognized as isoflavonoid (missing benzopyran fused benzene with appropriate aryl substituent)"

# Example usage (for testing):
if __name__ == '__main__':
    # Some test SMILES. (Note: many real isoflavonoids are large and decorated.)
    test_smiles = [
        "c1ccc2c(c1)oc(c2)c3ccccc3",  # simplified isoflavonoid (isoflavone type)
        "c1ccc2c(c1)oc(c2)C(c3ccccc3)",  # simplified isoflavonoid (isoflavan type)
        "O=C1C(=O)C2=CC=CC=C2O1",   # not isoflavonoid (coumarin)
    ]
    for smi in test_smiles:
        res, reason = is_isoflavonoid(smi)
        print("SMILES:", smi)
        print("Result:", res)
        print("Reason:", reason)
        print("="*60)