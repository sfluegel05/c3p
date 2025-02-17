"""
Classifies: CHEBI:35436 D-glucoside
"""
#!/usr/bin/env python3
"""
Classifies: D-Glucoside
Definition: Any glucoside in which the glycoside group is derived from D‐glucose.
For our purposes a positive is declared if we can find a six‐membered (pyranose) ring 
composed of five carbons and one oxygen that:
  (i) bears at least one exocyclic –CH2OH (or acyl‐modified CH2OH) branch, and 
 (ii) is linked via at least one external oxygen (the glycosidic bond) to a fragment 
     outside that ring.
Note: Many borderline cases exist.
"""
from rdkit import Chem

def is_D_glucoside(smiles: str):
    """
    Determines whether a molecule is a D-glucoside.
    
    Algorithm:
      1. Parse the SMILES string.
      2. Loop over all rings. For each ring of exactly 6 atoms:
           a. Confirm the ring contains exactly 1 oxygen and 5 carbons.
           b. Look for at least one ring carbon that has a substituent which is 
              a carbon (CH2 or longer) having at least one oxygen neighbor (i.e. a 
              relaxed –CH2OH check – even acylated groups are accepted).
           c. Look for at least one ring carbon that is linked via an external oxygen 
              (a potential glycosidic bond) to another heavy atom. In a typical hydroxyl 
              group (–OH) the oxygen is attached only to one heavy atom (the sugar). 
              A bridging (glycosidic) oxygen will connect to the sugar ring and at least 
              one additional heavy atom.
      3. If a candidate ring passes these tests, we declare the molecule to be a 
         D-glucoside.
    
    Args:
      smiles (str): The SMILES string of the molecule.
    
    Returns:
      bool: True if a candidate D-glucopyranose ring attached via a (likely) glycosidic bond is found, 
            False otherwise.
      str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        # Consider only 6-membered rings (pyranoses)
        if len(ring) != 6:
            continue

        atoms_in_ring = [mol.GetAtomWithIdx(i) for i in ring]
        # Count atoms: pyranose ring should have 5 carbons and 1 oxygen.
        n_oxygens = sum(1 for a in atoms_in_ring if a.GetAtomicNum() == 8)
        n_carbons  = sum(1 for a in atoms_in_ring if a.GetAtomicNum() == 6)
        if n_oxygens != 1 or n_carbons != 5:
            continue  # not a glucopyranose-like ring

        # Check for the characteristic exocyclic CH2OH substituent.
        found_ch2oh = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # only consider ring carbons
            if atom.GetAtomicNum() != 6:
                continue
            # Look among neighbors that lie outside the ring
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # We accept an exocyclic carbon even if it is modified (e.g. acylated)
                if nbr.GetAtomicNum() == 6:
                    # Check that this exocyclic carbon has at least one oxygen in its neighbors
                    if any(nbr2.GetAtomicNum() == 8 and nbr2.GetIdx() not in ring for nbr2 in nbr.GetNeighbors()):
                        found_ch2oh = True
                        break
            if found_ch2oh:
                break
        if not found_ch2oh:
            continue  # No characteristic CH2OH branch found

        # Check for a glycosidic bond: at least one ring carbon must connect via an external oxygen
        # to another heavy atom. (A simple hydroxyl would only have 1 heavy neighbor.)
        glycosidic_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # Count the heavy (non-hydrogen) neighbors of this oxygen
                    heavy_neighbors = [n for n in nbr.GetNeighbors() if n.GetAtomicNum() > 1]
                    # In a free hydroxyl group, the oxygen typically has only one heavy neighbor.
                    # In a glycosidic bond, the oxygen bridges the sugar ring and at least one other heavy atom.
                    if len(heavy_neighbors) >= 2:
                        glycosidic_found = True
                        break
            if glycosidic_found:
                break

        if glycosidic_found:
            return True, "Contains D-glucopyranose ring attached via a glycosidic bond"

    return False, "No D-glucopyranose ring attached via a glycosidic bond found"

# Example usage:
if __name__ == '__main__':
    # Test with a positive example: beta-D-glucosyl-N-eicosanoylsphingosine.
    test_smiles = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
    result, reason = is_D_glucoside(test_smiles)
    print(result, reason)