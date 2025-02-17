"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: Any 1-benzopyran with an aryl substituent at position 3 (isoflavonoid)

An isoflavonoid is defined as a benzopyran (chromene) system fused with a benzene ring,
with the heterocycle’s carbon adjacent to the oxygen (position 3) substituted by an aryl group.
In the previous attempt, strict rules (e.g. exactly one oxygen in the ring and any attached aromatic group)
led to many false negatives and false positives.
The improved algorithm proceeds as follows:
  1. Collect all six‐membered aromatic rings.
  2. Select “pyran candidates” that contain at least one oxygen.
  3. For each pyran candidate, look for a fused ring (sharing at least 2 atoms) that is a benzene ring (six‐membered, all carbons).
  4. Within the pyran candidate, identify the oxygen atom and then among its two ring neighbours choose the one that is not part of the fused intersection.
     This atom is our candidate “C-3”.
  5. Check that this candidate carbon has at least one neighbor (exterior to the pyran) that belongs to a disjoint six‐membered aromatic ring composed solely of carbons.
If these conditions are met we classify the molecule as an isoflavonoid.
Note: This algorithm is heuristic and may still mis‐classify border cases.
"""
from rdkit import Chem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    
    An isoflavonoid is defined as any 1-benzopyran with an aryl substituent at position 3.
    This improved implementation:
      1. Finds all 6-membered aromatic rings.
      2. Selects candidate rings that contain at least one oxygen (relaxing the “exactly one” rule).
      3. For each candidate fused system (a pyran candidate fused with a benzene ring composed solely of carbons):
            a. Identifies the oxygen atom in the pyran.
            b. Chooses the pyran neighbor (ideally corresponding to C-3) that is not involved in the fusion.
            c. Checks if that carbon bears an external substituent that is part of a separate benzene ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as an isoflavonoid, False otherwise.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string: parsing failed."
    
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in molecule."
    
    # Gather all six-membered aromatic rings (as sets of atom indices)
    all_six_arom = []
    for ring in ring_info:
        if len(ring) == 6:
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                all_six_arom.append(set(ring))
    if not all_six_arom:
        return False, "No six-membered aromatic rings found."
    
    # Define helper: does ring (set of indices) consist solely of carbons?
    def is_all_carbon(ring_set):
        return all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring_set)
    
    # Get pyran candidates: any six-membered ring with at least one oxygen.
    pyran_candidates = [ring for ring in all_six_arom if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring)]
    if not pyran_candidates:
        return False, "No six-membered aromatic ring with an oxygen (pyran candidate) found."
    
    # Now iterate over pyran candidates and try to find a fused benzene
    for pyran in pyran_candidates:
        # For each possible other ring (benzene ring must be all carbons)
        for other in all_six_arom:
            if pyran == other:
                continue
            if not is_all_carbon(other):
                continue
            # Look for fusion: they must share at least 2 atoms.
            fusion = pyran.intersection(other)
            if len(fusion) < 2:
                continue
            
            # Now we have a fused pair: pyran and benzene. Our fused system is:
            fused_atoms = pyran.union(other)
            
            # In the pyran candidate, try to find the oxygen atom.
            oxy_atoms = [idx for idx in pyran if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
            if not oxy_atoms:
                continue  # Should not happen but be cautious
            
            # For each oxygen in the pyran, check its neighbors that are within the ring.
            for oxy_idx in oxy_atoms:
                oxy_atom = mol.GetAtomWithIdx(oxy_idx)
                # Get neighbors inside the pyran (they are candidates for positions adjacent to oxygen)
                pyran_neighbors = [nbr.GetIdx() for nbr in oxy_atom.GetNeighbors() if nbr.GetIdx() in pyran]
                if not pyran_neighbors:
                    continue
                # We expect one of these to be involved in the fusion and one that is free (candidate C-3)
                # Filter out those that appear in the fusion with the benzene ring.
                candidate_c3s = [idx for idx in pyran_neighbors if idx not in fusion]
                # If none is found, try all neighbors (if both are fused then we might be in a borderline case)
                if not candidate_c3s:
                    candidate_c3s = pyran_neighbors
                # For each candidate C-3, check if it bears an exocyclic substituent that is part of another aromatic benzene ring.
                for c3_idx in candidate_c3s:
                    c3_atom = mol.GetAtomWithIdx(c3_idx)
                    # Look at neighbors outside the pyran candidate
                    for nbr in c3_atom.GetNeighbors():
                        nbr_idx = nbr.GetIdx()
                        if nbr_idx in pyran:
                            continue  # skip atoms inside the pyran core
                        # For a true aryl substituent, we expect a single bond from c3 to the attached group:
                        bond = mol.GetBondBetweenAtoms(c3_idx, nbr_idx)
                        if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                            continue
                        # Now check whether this neighbor belongs to a six-membered aromatic ring that is disjoint from the fused system.
                        for ring in all_six_arom:
                            if nbr_idx in ring and is_all_carbon(ring):
                                # We want the attached aromatic ring to be largely separate:
                                if ring.isdisjoint(fused_atoms):
                                    return True, ("Molecule contains a fused benzopyran (chromene) core (a six-membered ring with oxygen fused to a benzene) with an aryl substituent at the candidate C-3 position, consistent with an isoflavonoid scaffold.")
    return False, ("Scaffold not recognized as isoflavonoid. The molecule may be missing a fused benzopyran core with an appropriate external aryl substituent at the expected position.")

# Example usage (for testing purposes):
if __name__ == '__main__':
    test_smiles = [
        "c1ccc2c(c1)oc(c2)c3ccccc3",  # simplified isoflavonoid (isoflavone type)
        "c1ccc2c(c1)oc(c2)C(c3ccccc3)",  # simplified isoflavonoid (isoflavan type)
        "O=C1C(=O)C2=CC=CC=C2O1",       # not an isoflavonoid (coumarin)
    ]
    for smi in test_smiles:
        res, reason = is_isoflavonoid(smi)
        print("SMILES:", smi)
        print("Result:", res)
        print("Reason:", reason)
        print("="*60)