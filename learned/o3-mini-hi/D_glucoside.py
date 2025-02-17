"""
Classifies: CHEBI:35436 D-glucoside
"""
#!/usr/bin/env python3
"""
Classifies: D-Glucoside
Definition: Any glucoside in which the glycoside group is derived from D‐glucose.
For our purposes a positive is declared if we can find a six‐membered sugar ring
(i.e. a pyranose) composed of five carbons and one oxygen that (i) bears a “CH2OH”
branch (here accepted in a relaxed way so that even acylated –CH2OH groups count) and
(ii) is linked via at least one external oxygen (the glycosidic bond) to another fragment.
Note: Many borderline cases exist.
"""
from rdkit import Chem

def is_D_glucoside(smiles: str):
    """
    Determines whether a molecule is a D-glucoside.
    
    The algorithm:
      1. Parse the SMILES string.
      2. Loop over all rings. For each ring of exactly 6 atoms:
           a. Confirm the ring contains exactly 1 oxygen and 5 carbons.
           b. Look for at least one ring carbon that has a substituent which is a carbon (exocyclic)
              that in turn is bonded to at least one oxygen (to allow even modified –CH2OH groups).
          c. Look for at least one ring carbon that is connected (via a single oxygen)
             outside the ring to another fragment (the glycosidic bond).
      3. If such a candidate ring is found, we assume that it represents a D-glucopyranose unit.
    
    Args:
      smiles (str): The SMILES string of the molecule.
    
    Returns:
      bool: True if a candidate D-glucopyranose ring with a glycosidic bond is found, False otherwise.
      str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo()
    
    # Loop over all rings
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue  # only consider 6-membered rings (pyranoses)
        
        # Get atoms belonging to the ring; count carbons and oxygens.
        atoms_in_ring = [mol.GetAtomWithIdx(i) for i in ring]
        n_oxygens = sum(1 for a in atoms_in_ring if a.GetAtomicNum() == 8)
        n_carbons  = sum(1 for a in atoms_in_ring if a.GetAtomicNum() == 6)
        if n_oxygens != 1 or n_carbons != 5:
            continue  # not a pyranose-like ring
        
        # Check for the characteristic CH2OH branch.
        found_ch2oh = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            # Look at neighbors outside the ring
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 6:
                    # Relaxed test: allow exocyclic carbon even if substituted.
                    # If it has at least one oxygen neighbor (even if acylated) we treat it as a candidate.
                    if any(n.GetAtomicNum() == 8 for n in nbr.GetNeighbors() if n.GetIdx() not in ring):
                        found_ch2oh = True
                        break
            if found_ch2oh:
                break
        if not found_ch2oh:
            continue  # this ring does not display a proper CH2OH-like substituent
        
        # Check for a glycosidic bond: one of the ring carbons should be connected via an external oxygen.
        glycosidic_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # Check that this oxygen connects further to at least one atom not in the ring.
                    if any(n.GetIdx() not in ring and n.GetIdx() != atom.GetIdx() for n in nbr.GetNeighbors()):
                        glycosidic_found = True
                        break
            if glycosidic_found:
                break
        
        if glycosidic_found:
            return True, "Contains D-glucopyranose ring attached via a glycosidic bond"
    
    return False, "No D-glucopyranose ring attached via a glycosidic bond found"

# Example usage:
if __name__ == '__main__':
    # Test with one of the positive examples: beta-D-glucosyl-N-eicosanoylsphingosine.
    test_smiles = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
    result, reason = is_D_glucoside(test_smiles)
    print(result, reason)