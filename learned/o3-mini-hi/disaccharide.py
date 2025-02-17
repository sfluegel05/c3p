"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: Disaccharide
A disaccharide is defined as a compound in which two monosaccharides are joined by a glycosidic bond.
This algorithm works in two stages:
  1) Identify candidate monosaccharide rings: we search for 5- or 6-membered rings with exactly one ring oxygen.
  2) Look for a glycosidic linkage by searching for an oxygen atom whose neighbors include atoms found in two different candidate rings.
If both conditions are met, we return classification as disaccharide.
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    
    A disaccharide is assumed to contain exactly two monosaccharide rings (typically 5- or 6-membered rings
    with one ring oxygen) connected by a glycosidic linkage.
    
    The glycosidic linkage is identified by locating an oxygen atom that bridges an atom from one candidate sugar ring
    to an atom from the other candidate sugar ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a disaccharide, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # Each ring is represented as a tuple of atom indices.
    
    candidate_rings = []
    # Identify candidate sugar rings: 
    #   only rings of size 5 or 6 that have exactly one oxygen (heteroatom) in the ring.
    for ring in atom_rings:
        if len(ring) not in (5, 6):
            continue
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Count oxygen atoms that are actually part of the ring.
        ring_oxygens = [atom for atom in atoms_in_ring if atom.GetAtomicNum() == 8]
        if len(ring_oxygens) == 1:
            candidate_rings.append(set(ring))
    
    if len(candidate_rings) != 2:
        return False, f"Found {len(candidate_rings)} candidate sugar ring(s); exactly 2 are needed for a disaccharide."
    
    ring1, ring2 = candidate_rings[0], candidate_rings[1]
    
    # Look for a glycosidic linkage:
    # Instead of starting from a candidate carbon trying to reach a bridging oxygen,
    # we iterate over all oxygen atoms in the molecule. If any oxygen has at least one neighbor in ring1
    # and at least one (other) neighbor in ring2, then we accept that as a glycosidic linkage.
    glycosidic_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue  # only consider oxygen atoms
        neigh_rings = set()
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in ring1:
                neigh_rings.add("ring1")
            if nbr_idx in ring2:
                neigh_rings.add("ring2")
        if "ring1" in neigh_rings and "ring2" in neigh_rings:
            glycosidic_found = True
            break
    
    if glycosidic_found:
        return True, "Contains exactly two monosaccharide rings joined by a glycosidic bond."
    else:
        return False, "No glycosidic linkage bridging the two sugar rings was found."

# Example usage (for testing):
if __name__ == "__main__":
    # Provided test disaccharide: alpha-L-Fucp-(1->6)-alpha-D-Glcp
    test_smiles = "O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)CO[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C"
    result, reason = is_disaccharide(test_smiles)
    print(result, reason)