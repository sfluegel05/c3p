"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: Disaccharide
A disaccharide is defined as a compound in which two monosaccharides (typically in the form of a furanose or pyranose ring)
are joined by a glycosidic bond.
The algorithm works in two stages:
  1) Identify candidate monosaccharide rings by examining each ring (from rdkit's ring info) for being either a 5-membered (furanose)
     ring with exactly 1 oxygen and 4 carbons, or a 6-membered (pyranose) ring with exactly 1 oxygen and 5 carbons.
  2) Identify a glycosidic linkage by finding an oxygen (outside of the rings) that connects (via bonds) to at least one carbon in each of the two rings.
If both conditions are met then the molecule is classified as a disaccharide.
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    
    A disaccharide should contain exactly two candidate sugar rings (either 5-membered or 6-membered with one ring oxygen)
    connected by a glycosidic linkage (an oxygen bridging carbons from each ring).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a disaccharide, False otherwise.
        str: A reason for the classification result.
    """
    # Convert SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Retrieve ring information from the molecule
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # each ring is a tuple of atom indices

    candidate_rings = []  # to store sets of atom indices that look like sugar rings
    
    # Check each ring to see if it qualifies as a monosaccharide ring.
    for ring in rings:
        ring_atoms = list(ring)
        # Only consider rings of size 5 or 6
        if len(ring_atoms) not in (5, 6):
            continue
        # Count oxygen and carbon atoms in the ring.
        o_count = 0
        c_count = 0
        for idx in ring_atoms:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                o_count += 1
            elif atom.GetAtomicNum() == 6:
                c_count += 1
        # For a 5-membered ring (furanose): expect 1 oxygen and 4 carbons.
        if len(ring_atoms) == 5 and o_count == 1 and c_count == 4:
            candidate_rings.append(set(ring_atoms))
        # For a 6-membered ring (pyranose): expect 1 oxygen and 5 carbons.
        elif len(ring_atoms) == 6 and o_count == 1 and c_count == 5:
            candidate_rings.append(set(ring_atoms))
    
    # We expect exactly 2 candidate sugar rings for a disaccharide.
    if len(candidate_rings) != 2:
        return False, f"Found {len(candidate_rings)} candidate sugar ring(s); exactly 2 are needed for a disaccharide."
    
    ring1, ring2 = candidate_rings

    # Now look for the glycosidic bridge:
    # We search for an oxygen atom that is not part of either ring, and that is bonded to at least one carbon atom in each ring.
    glyco_found = False
    for atom in mol.GetAtoms():
        # Work only with oxygen atoms.
        if atom.GetAtomicNum() != 8:
            continue
        # Skip oxygens that are already part of one of the rings.
        if atom.GetIdx() in ring1 or atom.GetIdx() in ring2:
            continue
        
        # Check neighbors: we want at least one neighbor in ring1 and one in ring2 that is a carbon.
        neighbors = atom.GetNeighbors()
        in_ring1 = False
        in_ring2 = False
        for nbr in neighbors:
            if nbr.GetAtomicNum() != 6:   # only consider carbon neighbors
                continue
            if nbr.GetIdx() in ring1:
                in_ring1 = True
            if nbr.GetIdx() in ring2:
                in_ring2 = True
        if in_ring1 and in_ring2:
            glyco_found = True
            break
    
    if not glyco_found:
        return False, "No glycosidic linkage bridging the two sugar rings was found."
    
    return True, "Contains exactly two monosaccharide rings joined by a glycosidic bond."

# Example usage for testing (uncomment the code below to run tests):
if __name__ == "__main__":
    # Example: alpha-L-Fucp-(1->6)-alpha-D-Glcp
    test_smiles = "O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)CO[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C"
    result, reason = is_disaccharide(test_smiles)
    print(result, reason)