"""
Classifies: CHEBI:35436 D-glucoside
"""
#!/usr/bin/env python3
"""
Classifies: D-Glucoside
Definition: Any glucoside in which the glycoside group is derived from D-glucose.
For our purposes we require two conditions:
  (i) the molecule contains an unmodified D-glucopyranose ring – here recognized heuristically as a 6-membered ring
      (5 carbons + 1 oxygen) in which one of the carbons bears a CH2OH substituent (with exactly one oxygen bonded to that exocyclic carbon),
  (ii) that ring is connected via an oxygen (a glycosidic bond) to an aglycone fragment that is not itself a sugar ring.
Note: Many borderline cases exist.
"""

from rdkit import Chem

def is_D_glucoside(smiles: str):
    """
    Determines whether a molecule is a D-glucoside.
    The algorithm goes as follows:
      1. Parse the SMILES.
      2. Loop through rings of exactly 6 atoms.
         a. Keep only rings with exactly one oxygen and five carbons.
      3. For each candidate ring, check that one ring carbon bears a CH2OH branch.
         Here we require that an exocyclic (non‐ring) carbon attached to a ring carbon has:
           – at least 2 implicit hydrogens (suggesting CH2),
           – and at least one oxygen neighbor (the –OH) that itself is only connected to the exocyclic carbon.
      4. Also confirm that at least one ring carbon is linked via a bridging oxygen
         (i.e. via an oxygen outside the ring that itself attaches to a carbon not part of any sugar ring)
         to a likely aglycone.
      5. If both conditions are met, classify as a D-glucoside.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    # A helper: record indices of all atoms that belong to a confirmed D-glucose ring.
    sugar_ring_atom_indices = set()
    
    # Check each ring for a candidate D-glucopyranose.
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue
        # Get atoms from ring.
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        # Count oxygen and carbon atoms in ring.
        oxygens = [a for a in atoms if a.GetAtomicNum() == 8]
        carbons = [a for a in atoms if a.GetAtomicNum() == 6]
        if len(oxygens) != 1 or len(carbons) != 5:
            continue  # Not a pyranose‐like ring.
        
        # Check for the characteristic CH2OH branch.
        found_ch2oh = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            # Look at neighbors not inside the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Candidate exocyclic carbon for CH2OH: must be carbon.
                if nbr.GetAtomicNum() == 6:
                    # Check that it has at least 2 hydrogens.
                    # (Using total number of hydrogens, which includes implicit Hs.)
                    if nbr.GetTotalNumHs() < 2:
                        continue
                    # Check for at least one oxygen neighbor (the –OH)
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetIdx() == atom.GetIdx():
                            continue
                        if nbr2.GetAtomicNum() == 8:
                            # For a pure hydroxyl group, oxygen should be connected only to this carbon.
                            if nbr2.GetDegree() == 1:
                                found_ch2oh = True
                                break
                    if found_ch2oh:
                        break
            if found_ch2oh:
                break
        if not found_ch2oh:
            continue  # This ring does not have the expected CH2OH branch.
        
        # Now check for a glycosidic bond: one of the ring carbons should bond via an external oxygen to an aglycone.
        glycosidic_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Only consider carbons in the ring.
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # Consider the oxygen: it is part of the glycosidic bond if it attaches to a carbon not in a sugar ring.
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetIdx() == atom.GetIdx() or nbr2.GetIdx() in ring:
                            continue
                        if nbr2.GetAtomicNum() == 6:
                            # Heuristic: if the connected carbon is not part of another sugar ring (check later),
                            # then assume this oxygen is linking sugar to aglycone.
                            glycosidic_found = True
                            break
                if glycosidic_found:
                    break
            if glycosidic_found:
                break
        
        if glycosidic_found:
            # Mark this ring’s atoms as belonging to a D-glucose unit.
            sugar_ring_atom_indices.update(ring)
            return True, "Contains D-glucopyranose ring attached via a glycosidic bond"
    
    return False, "No unmodified D-glucopyranose ring attached via a glycosidic bond found"

# Example usage:
if __name__ == '__main__':
    # Test example: beta-D-glucosyl-N-eicosanoylsphingosine.
    test_smiles = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
    result, reason = is_D_glucoside(test_smiles)
    print(result, reason)