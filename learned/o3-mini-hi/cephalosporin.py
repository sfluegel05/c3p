"""
Classifies: CHEBI:23066 cephalosporin
"""
"""
Classifies: cephalosporin
Definition: A beta‐lactam antibiotic whose core comprises a fused bicyclic ring system in which a 
4‐membered beta‐lactam ring (with one ring‐bound nitrogen and a carbonyl group) is fused to a 6‐membered 
dihydrothiazine ring (with at least one sulfur). This program detects candidate rings and then verifies that 
a 4-membered and a 6-membered ring are fused (sharing exactly two atoms).
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    The function checks for the presence of a fused beta-lactam (4-membered) ring and
    dihydrothiazine (6-membered) ring by examining the rings in the molecule.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a cephalosporin core is detected, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse SMILES into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Ensure that ring information is computed
    # (GetRingInfo() is automatically computed but can be forced using GetSymmSSSR)
    rings = mol.GetRingInfo().AtomRings()
    
    four_membered_candidates = []
    six_membered_candidates = []
    
    # Loop over all rings in the molecule
    for ring in rings:
        ring_set = set(ring)
        if len(ring) == 4:
            # For a beta-lactam ring, we need exactly one ring nitrogen and at least one carbon with a carbonyl group.
            n_count = 0
            has_carbonyl = False
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 7:
                    n_count += 1
                if atom.GetAtomicNum() == 6:
                    # Check if any double bond from this carbon goes to an oxygen
                    for bond in atom.GetBonds():
                        if bond.GetBondType() == rdchem.BondType.DOUBLE:
                            nbr = bond.GetOtherAtom(atom)
                            if nbr.GetAtomicNum() == 8:
                                has_carbonyl = True
            if n_count == 1 and has_carbonyl:
                four_membered_candidates.append(ring_set)
                
        elif len(ring) == 6:
            # A dihydrothiazine ring should contain at least one sulfur atom.
            if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 16 for idx in ring):
                six_membered_candidates.append(ring_set)
    
    # Now check if any 4-membered candidate ring is fused with any 6-membered candidate ring (sharing exactly two atoms)
    for four_ring in four_membered_candidates:
        for six_ring in six_membered_candidates:
            common_atoms = four_ring.intersection(six_ring)
            if len(common_atoms) == 2:
                return True, "Cephalosporin core detected: fused 4-membered beta-lactam and 6-membered dihydrothiazine rings found."
    
    return False, "No cephalosporin core structure detected"

# Example usage/testing (uncomment the code below for testing purposes):
# if __name__ == '__main__':
#     test_smiles = "[H][C@]12SCC(C)=C(N1C(=O)[C@H]2N)C(O)=O"  # 7beta-aminodeacetoxycephalosporanic acid
#     result, reason = is_cephalosporin(test_smiles)
#     print("Is cephalosporin:", result)
#     print("Reason:", reason)