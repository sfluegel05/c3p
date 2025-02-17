"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: Disaccharide
A disaccharide is defined as a compound in which two monosaccharides are joined by a glycosidic bond.
This algorithm works in two stages:
  1) Identify candidate monosaccharide rings: we look for 5- or 6-membered rings that contain exactly one ring oxygen.
     (We allow some substituents to be present).
  2) Search for a glycosidic linkage by looking for an anomeric carbon (a carbon in one sugar ring having an oxygen neighbor outside the ring)
     that is connected through an oxygen atom to a carbon in the other ring.
If both conditions are met, we return classification as disaccharide.
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    
    A disaccharide is assumed to contain exactly two monosaccharide rings (typically 5- or 6-membered rings
    with one ring oxygen) connected by a glycosidic linkage. 
    The glycosidic linkage is identified by checking if a ring carbon (anomeric carbon candidate) in one ring
    is connected to an oxygen that in turn is bonded to a carbon in the other ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a disaccharide, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get fused ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # tuple of tuples: each tuple is atom indices of that ring
    
    candidate_rings = []
    # Identify candidate sugar rings: 5- or 6-membered rings that contain exactly one ring oxygen.
    for ring in atom_rings:
        if len(ring) not in (5, 6):
            continue  # only consider rings of size 5 or 6
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Count the number of oxygens that are members of the ring.
        ring_oxygens = [atom for atom in atoms_in_ring if atom.GetAtomicNum() == 8]
        # We expect one ring oxygen (the heterocycle) and the rest carbons (or possibly substituted with H, but not with another heteroatom)
        # Allow candidate if exactly one oxygen in the ring.
        if len(ring_oxygens) == 1:
            candidate_rings.append(set(ring))
    
    if len(candidate_rings) != 2:
        return False, f"Found {len(candidate_rings)} candidate sugar ring(s); exactly 2 are needed for a disaccharide."
    
    ring1, ring2 = candidate_rings[0], candidate_rings[1]
    
    # Look for a glycosidic linkage.
    # We search for an anomeric carbon in one ring. Many sugar rings have a carbon with a bond 
    # to an oxygen that is not part of the ring (the glycosidic oxygen).
    glycosidic_found = False
    for ringA, ringB in ((ring1, ring2), (ring2, ring1)):
        # For each candidate sugar ring (ringA), look at its carbon atoms as potential anomeric centers.
        for idx in ringA:
            atom = mol.GetAtomWithIdx(idx)
            # Focus on carbons only.
            if atom.GetAtomicNum() != 6:
                continue
            # Look through neighbors for an oxygen that is not (or not exclusively) part of ringA.
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() != 8:
                    continue
                nbr_idx = nbr.GetIdx()
                # Check: If this oxygen is not in ringA, or even if it is in ringA but has extra bonds, it can be bridging.
                # We then examine the neighbors of this oxygen to see if one of them belongs to ringB.
                # (Allow the possibility that the glycosidic oxygen is formally "in" a ring but still acts as a bridge.)
                # Do not consider the current carbon (atom) as the second neighbor.
                for onbr in nbr.GetNeighbors():
                    if onbr.GetIdx() == idx:
                        continue
                    if onbr.GetIdx() in ringB:
                        glycosidic_found = True
                        break
                if glycosidic_found:
                    break
            if glycosidic_found:
                break
        if glycosidic_found:
            break
    
    if glycosidic_found:
        return True, "Contains exactly two monosaccharide rings joined by a glycosidic bond."
    else:
        return False, "No glycosidic linkage found that connects the two sugar rings."

# Example usage (for testing):
if __name__ == "__main__":
    # Provided example disaccharide: alpha-L-Fucp-(1->6)-alpha-D-Glcp
    test_smiles = "O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)CO[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C"
    result, reason = is_disaccharide(test_smiles)
    print(result, reason)