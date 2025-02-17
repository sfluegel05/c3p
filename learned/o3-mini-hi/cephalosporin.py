"""
Classifies: CHEBI:23066 cephalosporin
"""
"""
Classifies: cephalosporin
Definition: A beta–lactam antibiotic having a fused bicyclic ring system where a 4–membered beta–lactam ring (with an amide group) is fused with a 6–membered dihydrothiazine ring (containing at least one sulfur atom).
"""

from rdkit import Chem

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    Cephalosporins have a fused bicyclic system consisting of a 4-membered beta-lactam ring (with a nitrogen and a carbonyl) 
    fused to a 6-membered dihydrothiazine ring (with at least one sulfur atom).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a cephalosporin, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information (list of tuples of atom indices for each ring)
    ring_info = mol.GetRingInfo().AtomRings()
    
    beta_lactam_rings = []  # 4-membered rings with N and a carbonyl C
    thiazine_rings = []     # 6-membered rings with at least one S
    
    # Helper: check if a given ring (list of atom indices) is a beta-lactam ring candidate.
    def is_beta_lactam_ring(ring):
        if len(ring) != 4:
            return False
        has_N = False
        has_carbonyl = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() == "N":
                has_N = True
            if atom.GetSymbol() == "C":
                # Check if carbon has a double bond to an oxygen (exocyclic carbonyl)
                for bond in atom.GetBonds():
                    # Get the neighbor atom that is not in the ring (or even if in ring, we consider exocyclic carbonyl bonds)
                    nbr = bond.GetOtherAtom(atom)
                    # Check for a double bond to oxygen
                    if nbr.GetSymbol() == "O" and bond.GetBondTypeAsDouble() == 2.0:
                        has_carbonyl = True
                        break
                # If one carbon already qualifies, no need to check further for carbonyl
                if has_carbonyl:
                    continue
        return has_N and has_carbonyl

    # Helper: check if a given ring is a dihydrothiazine ring candidate.
    def is_dihydrothiazine_ring(ring):
        if len(ring) != 6:
            return False
        # Look for at least one sulfur in the ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() == "S":
                return True
        return False

    # Collect beta-lactam and dihydrothiazine (6-membered with S) candidate rings.
    for ring in ring_info:
        if is_beta_lactam_ring(ring):
            beta_lactam_rings.append(set(ring))
        if is_dihydrothiazine_ring(ring):
            thiazine_rings.append(set(ring))
            
    if not beta_lactam_rings:
        return False, "No beta-lactam ring (4-membered ring with an amide group) found"
    if not thiazine_rings:
        return False, "No 6-membered ring with a sulfur atom (dihydrothiazine ring) found"
        
    # Check for fusion: a cephalosporin must have one beta-lactam ring fused with one dihydrothiazine ring
    for beta_ring in beta_lactam_rings:
        for thiazine_ring in thiazine_rings:
            # Fused rings should share exactly two atoms
            common_atoms = beta_ring.intersection(thiazine_ring)
            if len(common_atoms) == 2:
                return True, "Fused beta-lactam (4-membered) and dihydrothiazine (6-membered) ring system detected"
    
    return False, "No appropriate fused ring system (shared 2 atoms between a beta-lactam and a dihydrothiazine ring) found"

# Example usage (for testing purposes):
if __name__ == "__main__":
    # 7beta-aminodeacetoxycephalosporanic acid example SMILES
    test_smiles = "[H][C@]12SCC(C)=C(N1C(=O)[C@H]2N)C(O)=O"
    is_ceph, reason = is_cephalosporin(test_smiles)
    print("Is cephalosporin:", is_ceph)
    print("Reason:", reason)