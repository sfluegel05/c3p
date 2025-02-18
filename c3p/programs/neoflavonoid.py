"""
Classifies: CHEBI:71971 neoflavonoid
"""
"""
Classifies: Neoflavonoid
Definition: Any 1-benzopyran with an aryl substituent at position 4.
A neoflavonoid is defined here as a fused bicyclic system where
  - a 6-membered heterocycle (the chromene/coumarin core) containing exactly one oxygen is fused
    (i.e. sharing exactly 2 atoms, none of which is oxygen) with
  - an aromatic benzene ring (a 6-membered ring of carbons flagged as aromatic).
Furthermore, one of the non-fused carbons of the heterocycle (typically the one adjacent to the oxygen)
must have an external aryl substituent – here, defined as being connected to a separate 6-membered aromatic ring.
Note: Being “on the 4–position” is assumed to be the non-shared neighbor of the oxygen.
"""

from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid.
    
    Strategy:
      1. Parse the molecule and extract ring information.
      2. Identify candidate heterocycle rings that:
           - are 6-membered,
           - contain exactly one oxygen (assumed to be part of the pyran ring),
         and candidate benzene rings that:
           - are 6-membered,
           - all atoms are aromatic carbons.
      3. For each pair (heterocycle, benzene ring) that share exactly two atoms – and none of the shared atoms is oxygen –
         assume that the heterocycle is our chromene/coumarin core.
      4. In that heterocycle, find the oxygen atom. Among its neighbors in the heterocycle, at least one is expected to be “non-fused”
         (i.e. not in the shared set). This atom is our candidate “position 4.”
      5. Check that this candidate carbon has at least one neighbor (outside the heterocycle) that is aromatic and that
         belongs to a 6-membered aromatic ring (i.e. a phenyl substituent).
      6. If such a match is found, return True with an explanation; otherwise return False.
    
    Args:
       smiles (str): SMILES representation of the molecule.
    
    Returns:
       bool: True if it matches the neoflavonoid definition, False otherwise.
       str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    
    # Gather candidate heterocycle rings (6-membered, with exactly one O) and candidate benzene rings
    candidate_heterocycles = []
    candidate_benzenes = []
    for ring in all_rings:
        if len(ring) != 6:
            continue
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        # For benzene, all atoms must be aromatic carbons.
        if all(a.GetIsAromatic() and a.GetAtomicNum() == 6 for a in atoms):
            candidate_benzenes.append(set(ring))
        # For a candidate chromene/coumarin ring, require exactly one oxygen and at least one carbon.
        o_count = sum(1 for a in atoms if a.GetAtomicNum() == 8)
        c_count = sum(1 for a in atoms if a.GetAtomicNum() == 6)
        if o_count == 1 and c_count >= 1:
            candidate_heterocycles.append(set(ring))
    
    if not candidate_heterocycles:
        return False, "No 6-membered heterocycle with exactly one oxygen found (candidate chromene/coumarin core missing)"
    if not candidate_benzenes:
        return False, "No 6-membered fully aromatic benzene ring found"
    
    # Try to find a fused pair: heterocycle and benzene sharing exactly 2 atoms (and neither shared atom is oxygen)
    for hetero in candidate_heterocycles:
        for benzen in candidate_benzenes:
            shared = hetero.intersection(benzen)
            if len(shared) != 2:
                continue
            # Make sure none of the shared atoms is oxygen (the oxygen in a chromene is not fused to the benzene)
            if any(mol.GetAtomWithIdx(i).GetAtomicNum() == 8 for i in shared):
                continue

            # Now assume that hetero is the chromene core. Find the oxygen atom in the heterocycle.
            oxygen_idx = None
            for i in hetero:
                if mol.GetAtomWithIdx(i).GetAtomicNum() == 8:
                    oxygen_idx = i
                    break
            if oxygen_idx is None:
                continue  # unexpected, but skip if none found
            
            oxygen_atom = mol.GetAtomWithIdx(oxygen_idx)
            # Find neighbors of oxygen that lie in the heterocycle.
            hetero_neighbors = [nbr.GetIdx() for nbr in oxygen_atom.GetNeighbors() if nbr.GetIdx() in hetero]
            # We expect one of these neighbors to be non-fused (i.e. not in the shared set) – candidate for position 4.
            candidate_pos4s = [idx for idx in hetero_neighbors if idx not in shared]
            if not candidate_pos4s:
                continue  # no non-fused neighbor of oxygen found in heterocycle
            
            # For each candidate position (there might be more than one if ring symmetry is broken)
            for pos_idx in candidate_pos4s:
                pos_atom = mol.GetAtomWithIdx(pos_idx)
                # Skip if this atom is not carbon.
                if pos_atom.GetAtomicNum() != 6:
                    continue
                # Check substituents attached to pos_atom that are NOT part of the heterocycle.
                for nbr in pos_atom.GetNeighbors():
                    if nbr.GetIdx() in hetero:
                        continue  # skip atoms belonging to the core
                    # We want an aryl substituent: first the neighbor should be aromatic.
                    if not nbr.GetIsAromatic():
                        continue
                    # Now check if this neighbor belongs to a 6-membered aromatic ring (phenyl ring)
                    for ring in all_rings:
                        if nbr.GetIdx() in ring and len(ring) == 6:
                            # Verify all atoms in that ring are aromatic carbons.
                            if all(mol.GetAtomWithIdx(j).GetIsAromatic() and mol.GetAtomWithIdx(j).GetAtomicNum() == 6 for j in ring):
                                return True, ("Found candidate 1-benzopyran core (heterocycle fused with benzene) "
                                              "with an aryl substituent at a candidate non-fused position "
                                              "(likely the 4-position)")
    return False, "No suitable candidate 1-benzopyran core with an external aryl substituent at the 4-position was identified"


# Example usage for testing:
if __name__ == "__main__":
    # Test with one of the provided SMILES strings
    test_smiles = "ClC1=C(O)C(=C(C(=O)OC)C=C1OC)C(=O)C2=C(O)C(C=3C4=C(OC(C3)(C)C)C(=CC(=C4)O)CC=C(C)C)=C(C)C=C2O"  # Pestaloficiol U
    result, reason = is_neoflavonoid(test_smiles)
    print("Is neoflavonoid?", result)
    print("Reason:", reason)