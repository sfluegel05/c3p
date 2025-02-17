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
    Algorithm improvements:
      - Use RDKit’s GetSymmSSSR instead of raw ring info, to limit to smallest rings.
      - For the beta-lactam ring (4-membered): require exactly one N and three C atoms (no others) and at least one of those carbons must have a double-bond to oxygen (exocyclic or in-ring).
      - For the dihydro(thia/oxa)zine ring (6-membered): require it to contain at least one S or O and not be fully aromatic (i.e. at least one atom is non aromatic).
      - The two rings are fused if they share exactly 2 atoms that are adjacent in each ring (i.e. the shared bond is contiguous).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a cephalosporin, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get the smallest rings using GetSymmSSSR
    ssr = Chem.GetSymmSSSR(mol)
    rings = [tuple(r) for r in ssr]
    
    # Helper: check if two indices are adjacent in a given ring (cyclic order)
    def are_adjacent_in_ring(ring, a, b):
        n = len(ring)
        for i in range(n):
            if (ring[i] == a and ring[(i+1)%n] == b) or (ring[i] == b and ring[(i+1)%n] == a):
                return True
        return False

    # Check for beta-lactam candidate: exactly 4 atoms, one nitrogen and three carbons.
    # Also ensure at least one carbon in the ring has a double bond (in or exocyclic) to an oxygen.
    def is_beta_lactam_ring(ring):
        if len(ring) != 4:
            return False
        n_count = 0
        carbon_indices = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            sym = atom.GetSymbol()
            if sym == "N":
                n_count += 1
            elif sym == "C":
                carbon_indices.append(idx)
            else:
                # if we have any non C or N atoms then it is unexpected.
                return False
        if n_count != 1 or len(carbon_indices) != 3:
            return False
        
        # Must find at least one carbon that is bound (by a double bond) to oxygen.
        carbonyl_found = False
        for idx in carbon_indices:
            atom = mol.GetAtomWithIdx(idx)
            for bond in atom.GetBonds():
                # bond double type
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    nbr = bond.GetOtherAtom(atom)
                    if nbr.GetSymbol() == "O":
                        carbonyl_found = True
                        break
            if carbonyl_found:
                break
        return carbonyl_found

    # Check for candidate 6-membered dihydro(thia/oxa)zine ring:
    # Must be exactly 6 atoms, contain at least one S or O, and not be fully aromatic.
    def is_dihydro_thia_oxa_zine_ring(ring):
        if len(ring) != 6:
            return False
        hetero_found = False
        non_aromatic_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() in ["S", "O"]:
                hetero_found = True
            if not atom.GetIsAromatic():
                non_aromatic_found = True
        return hetero_found and non_aromatic_found

    beta_lactam_rings = []
    thia_oxa_rings = []
    for ring in rings:
        if is_beta_lactam_ring(ring):
            beta_lactam_rings.append(ring)
        if is_dihydro_thia_oxa_zine_ring(ring):
            thia_oxa_rings.append(ring)
    
    if not beta_lactam_rings:
        return False, "No candidate beta-lactam ring (4-membered with 1N and 3C, and a carbonyl) found"
    if not thia_oxa_rings:
        return False, "No candidate 6-membered dihydro(thia/oxa)zine ring (with S or O and non-aromatic character) found"
    
    # For every candidate beta-lactam and thia/oxa ring, verify they fuse properly.
    # Proper fusion: share exactly 2 atoms and the shared atoms must be contiguous (i.e. adjacent) in both rings.
    for beta_ring in beta_lactam_rings:
        beta_set = set(beta_ring)
        for thia_ring in thia_oxa_rings:
            thia_set = set(thia_ring)
            common_atoms = beta_set.intersection(thia_set)
            if len(common_atoms) == 2:
                common_list = list(common_atoms)
                # Check that the common atoms are directly bonded
                bond = mol.GetBondBetweenAtoms(common_list[0], common_list[1])
                if bond is None:
                    continue
                if are_adjacent_in_ring(beta_ring, common_list[0], common_list[1]) and \
                   are_adjacent_in_ring(thia_ring, common_list[0], common_list[1]):
                    return True, "Fused beta-lactam (4-membered) and dihydro(thia/oxa)zine (6-membered) ring system detected"
    
    return False, "No fused ring system with proper beta-lactam and dihydro(thia/oxa)zine rings found"

# Example usage (for testing purposes only)
if __name__ == '__main__':
    # Test with 7beta-aminodeacetoxycephalosporanic acid SMILES
    test_smiles = "[H][C@]12SCC(C)=C(N1C(=O)[C@H]2N)C(O)=O"
    is_ceph, reason = is_cephalosporin(test_smiles)
    print("Is cephalosporin:", is_ceph)
    print("Reason:", reason)