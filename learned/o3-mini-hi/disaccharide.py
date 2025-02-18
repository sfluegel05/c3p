"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: Disaccharide
A disaccharide is defined as a compound in which two monosaccharides are joined by a glycosidic bond.
This algorithm makes two key checks:
  1) Identify candidate monosaccharide rings – we assume these are 5-membered (furanose) or 6-membered (pyranose)
     rings that contain exactly one ring oxygen and the remaining ring atoms are carbons.
  2) Check that exactly two such rings are found and that there is at least one exocyclic oxygen 
     that bonds to one atom in one ring and one atom in the other (i.e. a glycosidic linkage).
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.

    A disaccharide should be composed of exactly two monosaccharide rings (typically 5- or 6-membered rings
    with one ring oxygen and the rest carbons) connected by a glycosidic linkage, which is assumed here to be
    an exocyclic oxygen that is simultaneously bonded to atoms in both rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string into molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Extract ring information from the molecule
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # tuple of tuples; each tuple contains atom indices for that ring

    candidate_rings = []  # list of sets; each set is the indices of a candidate monosaccharide ring
    # Look for rings that are either 5- or 6-membered and contain exactly one oxygen with the rest carbons.
    for ring in atom_rings:
        if len(ring) not in (5, 6):
            continue  # only consider 5- or 6-membered rings
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        oxygen_count = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 8)
        carbon_count = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 6)
        # In a typical sugar ring, there is exactly one ring oxygen and the other atoms are carbons.
        if oxygen_count == 1 and carbon_count == (len(ring) - 1):
            candidate_rings.append(set(ring))

    # For a disaccharide we require exactly two candidate sugar rings.
    if len(candidate_rings) != 2:
        return False, f"Found {len(candidate_rings)} candidate sugar ring(s); exactly 2 are needed for a disaccharide."

    ring1, ring2 = candidate_rings[0], candidate_rings[1]

    # Look for a glycosidic linkage.
    # We assume that the glycosidic bond is an exocyclic oxygen (i.e. not in either ring)
    # that bridges atoms from ring1 and ring2.
    glycosidic_found = False
    for atom in mol.GetAtoms():
        # Consider only oxygen atoms
        if atom.GetAtomicNum() != 8:
            continue
        idx = atom.GetIdx()
        # Skip oxygen atoms that are part of one of the candidate rings;
        # the bridging oxygen should lie outside the ring.
        if idx in ring1 or idx in ring2:
            continue
        neighbors = atom.GetNeighbors()
        if len(neighbors) < 2:
            continue  # unlikely to be a bridging oxygen
        # Check if one neighbor belongs to ring1 and a different neighbor belongs to ring2.
        in_ring1 = any(nei.GetIdx() in ring1 for nei in neighbors)
        in_ring2 = any(nei.GetIdx() in ring2 for nei in neighbors)
        if in_ring1 and in_ring2:
            glycosidic_found = True
            break

    if glycosidic_found:
        return True, "Contains exactly two monosaccharide rings joined by a glycosidic bond."
    else:
        return False, "No glycosidic linkage found that connects the two sugar rings."

# Example usage (for testing):
if __name__ == "__main__":
    # Test with one of the known disaccharide SMILES strings
    test_smiles = "O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)CO[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C"
    result, reason = is_disaccharide(test_smiles)
    print(result, reason)