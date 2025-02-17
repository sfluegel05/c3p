"""
Classifies: CHEBI:23066 cephalosporin
"""
"""
Classifies: cephalosporin
Definition: A beta‐lactam antibiotic whose core comprises a fused bicyclic ring system in which a 4‐membered beta‐lactam 
ring (with a ring-bound nitrogen and a carbonyl-bearing carbon) is fused (sharing exactly two adjacent atoms) to a 
6‐membered (dihydro)thia/oxa‐zine ring. In most cases the 6‐membered ring will contain a sulfur atom; however, if a ring 
contains an oxygen instead then it may be an oxacephalosporin.
"""

from rdkit import Chem

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    The algorithm:
      (1) Parses the molecule.
      (2) Finds candidate beta-lactam rings: exactly 4-membered rings that have exactly one nitrogen 
          and three carbon atoms, one of which must have a double bond O (carbonyl) attached (exocyclic or in‐ring).
      (3) Finds candidate dihydro(thia/oxa)zine rings: exactly 6-membered rings that contain at least one sulfur or 
          oxygen and are not fully aromatic.
      (4) For each candidate pair, requires that the rings are fused (they share exactly 2 atoms) and that the two 
          common atoms are adjacent in the ring order (i.e. the shared bond is contiguous in both rings).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a cephalosporin, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # helper: check if two indices are adjacent in a ring (cyclic order)
    def are_adjacent_in_ring(ring, a, b):
        # ring is a tuple of atom indices (in order)
        n = len(ring)
        for i in range(n):
            if (ring[i] == a and ring[(i+1)%n] == b) or (ring[i] == b and ring[(i+1)%n] == a):
                return True
        return False

    # candidate beta-lactam ring: 4-membered ring with exactly one N and three C atoms,
    # and at least one of the carbons should have a double-bonded oxygen.
    def is_beta_lactam_ring(ring):
        if len(ring) != 4:
            return False
        n_count = 0
        c_atoms = []
        carbonyl_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            sym = atom.GetSymbol()
            if sym == "N":
                n_count += 1
            elif sym == "C":
                c_atoms.append(atom)
            else:
                # any other atom in a 4-membered beta-lactam would be unexpected
                return False
        if n_count != 1 or len(c_atoms) != 3:
            return False
        # Check that at least one carbon in this ring has a double bond to oxygen
        for c in c_atoms:
            for bond in c.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    nbr = bond.GetOtherAtom(c)
                    if nbr.GetSymbol() == "O":
                        # We allow the carbonyl oxygen to be exocyclic or even if the oxygen is in the ring.
                        carbonyl_found = True
                        break
            if carbonyl_found:
                break
        return carbonyl_found

    # candidate dihydro(thia/oxa)zine ring: exactly 6-membered ring with at least one heteroatom (S or O)
    # and not fully aromatic (i.e. “dihydro”)
    def is_dihydro_thia_oxa_zine_ring(ring):
        if len(ring) != 6:
            return False
        hetero_found = False
        all_aromatic = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() in ["S", "O"]:
                hetero_found = True
            if not atom.GetIsAromatic():
                all_aromatic = False
        # We require at least one S or O; however if the entire ring is aromatic it is unlikely to be a dihydro ring.
        return hetero_found and (not all_aromatic)

    ring_info = mol.GetRingInfo().AtomRings()
    beta_lactam_rings = []  # list of rings (as tuples)
    thia_oxa_rings = []     # list of rings (as tuples)

    for ring in ring_info:
        if is_beta_lactam_ring(ring):
            beta_lactam_rings.append(ring)
        if is_dihydro_thia_oxa_zine_ring(ring):
            thia_oxa_rings.append(ring)
            
    if not beta_lactam_rings:
        return False, "No candidate beta-lactam ring (4-membered with 1N and carbonyl-bearing C) found"
    if not thia_oxa_rings:
        return False, "No candidate 6-membered dihydro(thia/oxa)zine ring (with heteroatom S/O) found"
        
    # Check for fusion: a fused system must have exactly two common atoms that are also adjacent in both rings.
    for beta_ring in beta_lactam_rings:
        beta_set = set(beta_ring)
        for thia_ring in thia_oxa_rings:
            thia_set = set(thia_ring)
            common_atoms = beta_set.intersection(thia_set)
            if len(common_atoms) == 2:
                common_list = list(common_atoms)
                # check that the two common atoms are connected by a bond in the molecule
                bond = mol.GetBondBetweenAtoms(common_list[0], common_list[1])
                if bond is None:
                    continue
                # also verify that in each ring the common atoms are adjacent (taking ring cyclic order into account)
                if are_adjacent_in_ring(beta_ring, common_list[0], common_list[1]) and are_adjacent_in_ring(thia_ring, common_list[0], common_list[1]):
                    return True, "Fused beta-lactam (4-membered) and dihydro(thia/oxa)zine (6-membered) ring system detected"
    
    return False, "No fused ring system with a proper beta-lactam and dihydro(thia/oxa)zine ring (sharing 2 adjacent atoms) found"

# Example usage (for testing purposes):
if __name__ == '__main__':
    # example: 7beta-aminodeacetoxycephalosporanic acid
    test_smiles = "[H][C@]12SCC(C)=C(N1C(=O)[C@H]2N)C(O)=O"
    is_ceph, reason = is_cephalosporin(test_smiles)
    print("Is cephalosporin:", is_ceph)
    print("Reason:", reason)