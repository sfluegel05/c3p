"""
Classifies: CHEBI:23066 cephalosporin
"""
"""
Classifies: cephalosporin
Definition: A beta‐lactam antibiotic whose core comprises a fused bicyclic ring system in which a 4‐membered beta‐lactam 
ring (with one ring-bound nitrogen and three carbons – one of which bears a carbonyl group) is fused (sharing exactly two adjacent atoms) 
to a 6‐membered dihydro(thia/oxa)zine ring (which must contain at least one sulfur or oxygen and be non‐aromatic). The fusion must include 
the carbonyl carbon from the beta-lactam.
"""

from rdkit import Chem

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    Improvements over the previous attempt:
      - When examining candidate beta–lactam rings (4-membered with 1N and 3C), we record which carbon has a double bond to oxygen.
      - Later, when checking fusion with a candidate 6-membered dihydro(thia/oxa)zine ring, we require
        that the two shared atoms between the rings include this beta–lactam carbonyl (i.e. the carbon with the double-bonded O).
      - We retain the requirement that the 6–membered ring is non‐aromatic and contains at least one S or O.
      - The two rings are considered fused if they share exactly 2 atoms that are adjacent (i.e. forming a bond) in both rings.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a cephalosporin, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get the smallest rings using RDKit's GetSymmSSSR
    ssr = Chem.GetSymmSSSR(mol)
    rings = [tuple(r) for r in ssr]
    
    # Helper: check if two indices are adjacent in a given cyclic ring
    def are_adjacent_in_ring(ring, a, b):
        n = len(ring)
        for i in range(n):
            if (ring[i] == a and ring[(i+1)%n] == b) or (ring[i] == b and ring[(i+1)%n] == a):
                return True
        return False

    # Check candidate beta-lactam ring.
    # Returns a tuple (True, carbonyl_idx) if the ring qualifies or (False, None) if not.
    def check_beta_lactam_ring(ring):
        if len(ring) != 4:
            return (False, None)
        n_count = 0
        carbon_indices = []
        carbonyl_index = None
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            sym = atom.GetSymbol()
            if sym == "N":
                n_count += 1
            elif sym == "C":
                carbon_indices.append(idx)
            else:
                # Any atom other than C or N disqualifies the ring.
                return (False, None)
        if n_count != 1 or len(carbon_indices) != 3:
            return (False, None)
        
        # Look for the carbon with a double bond to an oxygen (could be exocyclic or in-ring)
        for idx in carbon_indices:
            atom = mol.GetAtomWithIdx(idx)
            for bond in atom.GetBonds():
                # Consider a double bond to oxygen
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    nbr = bond.GetOtherAtom(atom)
                    if nbr.GetSymbol() == "O":
                        carbonyl_index = idx
                        break
            if carbonyl_index is not None:
                break
        if carbonyl_index is None:
            return (False, None)
        return (True, carbonyl_index)

    # Check candidate 6-membered dihydro(thia/oxa)zine ring.
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

    # Identify candidate rings.
    beta_lactam_candidates = []  # list of tuples: (ring, carbonyl_index)
    thia_oxa_candidates = []      # list of rings (tuple of indices)
    
    for ring in rings:
        valid_beta, carbonyl_idx = check_beta_lactam_ring(ring)
        if valid_beta:
            beta_lactam_candidates.append((ring, carbonyl_idx))
        if is_dihydro_thia_oxa_zine_ring(ring):
            thia_oxa_candidates.append(ring)
    
    if not beta_lactam_candidates:
        return False, "No candidate beta-lactam ring (4-membered with 1N, 3C, and a carbonyl) found"
    if not thia_oxa_candidates:
        return False, "No candidate 6-membered dihydro(thia/oxa)zine ring (with S or O and non-aromatic) found"
    
    # Check for correct fusion.
    # The rings must share exactly 2 atoms that are adjacent in each ring,
    # and the beta-lactam carbonyl carbon must be among these shared atoms.
    for beta_ring, carbonyl_idx in beta_lactam_candidates:
        beta_set = set(beta_ring)
        for thia_ring in thia_oxa_candidates:
            thia_set = set(thia_ring)
            common_atoms = beta_set.intersection(thia_set)
            if len(common_atoms) == 2:
                # The shared atoms must be adjacent in both rings.
                common_list = list(common_atoms)
                bond = mol.GetBondBetweenAtoms(common_list[0], common_list[1])
                if bond is None:
                    continue
                if not (are_adjacent_in_ring(beta_ring, common_list[0], common_list[1]) and
                        are_adjacent_in_ring(thia_ring, common_list[0], common_list[1])):
                    continue
                # Now require that the beta-lactam carbonyl carbon is among the shared atoms.
                if carbonyl_idx in common_atoms:
                    return True, "Fused beta-lactam (4-membered) and dihydro(thia/oxa)zine (6-membered) ring system detected"
    
    return False, "No fused ring system with proper beta-lactam (with carbonyl) and dihydro(thia/oxa)zine rings found"

# Example usage (for testing purposes only)
if __name__ == '__main__':
    # Test with 7beta-aminodeacetoxycephalosporanic acid SMILES
    test_smiles = "[H][C@]12SCC(C)=C(N1C(=O)[C@H]2N)C(O)=O"
    is_ceph, reason = is_cephalosporin(test_smiles)
    print("Is cephalosporin:", is_ceph)
    print("Reason:", reason)