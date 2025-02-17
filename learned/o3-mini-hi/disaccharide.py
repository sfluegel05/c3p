"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: Disaccharide
A disaccharide is defined as a compound in which two monosaccharides are joined by a glycosidic bond.
The idea here is to:
  1) Identify monosaccharide rings – typically 5-membered (furanose) or 6-membered (pyranose)
     rings that contain exactly one oxygen atom in the ring.
  2) Detect a glycosidic linkage – an oxygen atom bridging a carbon of one candidate sugar ring to a carbon of another.
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.

    A disaccharide should have two monosaccharide rings (typically 5- or 6-membered rings with one ring oxygen)
    and a glycosidic bond connecting the two (an oxygen from one ring connected to a carbon of the other).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # tuple of tuples with atom indices

    # Identify candidate monosaccharide rings:
    # We assume sugar rings are 5-membered (furanose) or 6-membered (pyranose)
    # and contain exactly one oxygen atom (atomic number 8) in the ring.
    candidate_sugar_rings = []  # each candidate is a set of atom indices
    for ring in atom_rings:
        if len(ring) not in (5, 6):
            continue  # skip rings not size 5 or 6
        # Count oxygen atoms in the ring
        oxygen_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:  # oxygen
                oxygen_count += 1
        if oxygen_count == 1:
            candidate_sugar_rings.append(set(ring))
            
    if len(candidate_sugar_rings) < 2:
        return False, f"Found only {len(candidate_sugar_rings)} candidate sugar unit(s), need at least 2"

    # Now, search for a glycosidic linkage between two candidate sugar rings.
    # We'll look for an exocyclic oxygen on a candidate sugar ring that connects to a carbon in another candidate ring.
    glycosidic_found = False
    # We check every candidate sugar ring
    for i, ring1 in enumerate(candidate_sugar_rings):
        for idx in ring1:
            atom = mol.GetAtomWithIdx(idx)
            # We are interested in carbons in the sugar ring that can be anomeric carbons.
            if atom.GetAtomicNum() != 6:
                continue
            # Loop over neighbors of the carbon to see if there is an exocyclic oxygen (i.e., oxygen not in the ring)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() != 8:
                    continue
                # If the oxygen atom is also in the ring, it is the ring oxygen. We need an exocyclic oxygen.
                if neighbor.GetIdx() in ring1:
                    continue
                # Now check the neighbors of this exocyclic oxygen.
                # If one of its neighbors (other than our original carbon) belongs to a different candidate sugar ring, we have a glycosidic bond.
                for oxy_nb in neighbor.GetNeighbors():
                    if oxy_nb.GetIdx() == idx:
                        continue  # skip the current carbon
                    # Check if this neighboring atom belongs to a different candidate sugar ring.
                    for j, ring2 in enumerate(candidate_sugar_rings):
                        if j == i:
                            continue
                        if oxy_nb.GetIdx() in ring2:
                            glycosidic_found = True
                            break
                    if glycosidic_found:
                        break
                if glycosidic_found:
                    break
            if glycosidic_found:
                break
        if glycosidic_found:
            break

    if glycosidic_found:
        return True, "Contains two monosaccharide rings joined by a glycosidic bond."
    else:
        return False, "No glycosidic linkage found between two monosaccharide rings."

# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = "O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)CO[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C"
    result, reason = is_disaccharide(test_smiles)
    print(result, reason)