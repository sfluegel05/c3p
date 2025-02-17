"""
Classifies: CHEBI:35436 D-glucoside
"""
"""
Classifies: D-Glucoside
Definition: Any glucoside in which the glycoside group is derived from D-glucose.
For our purposes we require two conditions:
  (i) the molecule contains an unmodified D-glucopyranose ring – here recognized heuristically as a 6-membered ring
      (5 carbons + 1 oxygen) in which one of the carbons bears a -CH2OH substituent (typical for D-glucose) and 
  (ii) that that ring is attached to an external moiety via an oxygen (i.e. a glycosidic bond).
  
Note: In practice there are many ways to depict a sugar. Rather than write two extremely specific SMARTS patterns,
this program examines the ring(s) in the molecule for the simplest characteristics of a D-glucopyranose.
Some true D-glucosides may still be missed and some non-D-glucosides misclassified.
"""

from rdkit import Chem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    
    The algorithm proceeds as follows:
      1. Parse the SMILES.
      2. Look at each ring of exactly 6 atoms.
      3. For each 6-atom ring, check that it contains exactly one oxygen and five carbons.
      4. Within that ring, check that at least one carbon carries a CH2OH branch (a signature of D-glucose).
         (Here we require that one neighboring atom outside of the ring is a carbon that itself is
          attached to an oxygen as its only heavy-atom substituent.)
      5. Also check that one of the ring carbons is connected via an oxygen (external to the ring)
         to a carbon (i.e. the sugar is “glycosylated” rather than a free sugar).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a D-glucoside, otherwise False.
        str: A short explanation.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    # Loop over all rings in the molecule.
    for ring in ring_info.AtomRings():
        # We are only interested in 6-membered rings.
        if len(ring) != 6:
            continue
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        # Count ring atoms by element.
        o_in_ring = [a for a in atoms if a.GetAtomicNum() == 8]
        c_in_ring = [a for a in atoms if a.GetAtomicNum() == 6]
        if len(o_in_ring) != 1 or len(c_in_ring) != 5:
            continue  # not a pyranose type ring
        # Check for the characteristic CH2OH branch.
        # In D-glucopyranose the “CH2OH” is attached to one of the ring carbons.
        found_ch2oh = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            # Look for a neighbor (outside the ring) that is carbon.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring: 
                    continue
                if nbr.GetAtomicNum() == 6:
                    # Check that this exocyclic carbon (the candidate CH2) is connected predominantly to an oxygen.
                    # We require that aside from the bond to the in-ring carbon, its heavy atoms count is 1 (an oxygen).
                    heavy_nbrs = [x for x in nbr.GetNeighbors() if x.GetAtomicNum() > 1 and x.GetIdx() != atom.GetIdx()]
                    if len(heavy_nbrs) == 1 and heavy_nbrs[0].GetAtomicNum() == 8:
                        found_ch2oh = True
                        break
            if found_ch2oh:
                break
        if not found_ch2oh:
            continue  # This ring does not look like an unmodified D-glucopyranose.
        
        # Check that the sugar is attached to an aglycone via a glycosidic bond.
        # Look for a bond from one of the ring carbons to an oxygen that in turn is attached to a carbon outside the ring.
        glyco_link = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # Check if this oxygen connects further to a carbon outside the ring.
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetIdx() != atom.GetIdx() and nbr2.GetIdx() not in ring and nbr2.GetAtomicNum() == 6:
                            glyco_link = True
                            break
                if glyco_link:
                    break
            if glyco_link:
                break
        
        if glyco_link:
            return True, "Contains D-glucopyranose ring attached via a glycosidic bond"
    return False, "No unmodified D-glucopyranose ring attached via a glycosidic bond found"


# Example usage
if __name__ == '__main__':
    # Test on beta-D-glucosyl-N-eicosanoylsphingosine (one known true example)
    test_smiles = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
    result, reason = is_D_glucoside(test_smiles)
    print(result, reason)