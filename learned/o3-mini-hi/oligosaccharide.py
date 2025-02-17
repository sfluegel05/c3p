"""
Classifies: CHEBI:50699 oligosaccharide
"""
#!/usr/bin/env python
"""
Classifies: Oligosaccharide
An oligosaccharide is defined as a compound in which at least two monosaccharide units (ring systems commonly 
found in sugars) are joined by glycosidic (typically ether) linkages.
This heuristic implementation identifies candidate sugar rings (of size 5 or 6 with one ring oxygen) 
and then looks for oxygen atoms that bridge atoms from two different rings.
"""

from rdkit import Chem

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string. It uses simple heuristics:
      1. Identify candidate sugar rings: cyclic substructures of size 5 (furanose) or 6 (pyranose)
         that contain exactly one oxygen atom and no aromatic atoms.
      2. Search for glycosidic oxygen linkages – non-ring oxygen atoms having two neighbors that belong 
         to different candidate sugar rings.
      3. If there are at least two sugar rings and at least one glycosidic linkage bridging them,
         then the molecule is classified as an oligosaccharide.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an oligosaccharide, else False.
        str: A reason explaining the outcome.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # tuple of tuples of atom indices
    
    # List to hold candidate sugar rings. Each entry is a set of atom indices.
    sugar_rings = []
    
    # For mapping: atom idx -> list of sugar ring indices that include this atom.
    atom_to_sugar_ring = {}
    
    # Iterate over each ring found in the molecule.
    for ring in atom_rings:
        ring_size = len(ring)
        # We are interested in rings of size 5 (furanose) or 6 (pyranose)
        if ring_size not in (5, 6):
            continue
        
        # Count the number of oxygen atoms in this ring and also check that none of the atoms are aromatic.
        oxygen_count = 0
        aromatic_flag = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                oxygen_count += 1
            if atom.GetIsAromatic():
                aromatic_flag = True
                break

        # In a typical carbohydrate ring, there is exactly one oxygen atom.
        if aromatic_flag:
            continue
        if oxygen_count != 1:
            continue

        # We can optionally try to demand that several substituents (hydroxyl groups) are attached,
        # but for this heuristic, we count the ring as a candidate sugar ring.
        ring_set = set(ring)
        sugar_rings.append(ring_set)
        ring_index = len(sugar_rings) - 1
        for idx in ring_set:
            atom_to_sugar_ring.setdefault(idx, []).append(ring_index)
    
    if len(sugar_rings) < 2:
        return False, f"Found {len(sugar_rings)} candidate sugar ring(s); need at least 2 for oligosaccharide"
    
    # Now look for glycosidic linkages.
    # The glycosidic oxygen should not be part of the candidate sugar ring itself.
    glycosidic_linkage_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        # Skip if the oxygen is part of any candidate sugar ring (it might be the ring oxygen)
        if atom.GetIdx() in atom_to_sugar_ring:
            continue
        
        # Consider only oxygens that have exactly two neighbors
        if atom.GetDegree() != 2:
            continue
        
        neighbor_sugar_rings = set()
        for neigh in atom.GetNeighbors():
            n_idx = neigh.GetIdx()
            if n_idx in atom_to_sugar_ring:
                # A neighbor may belong to more than one sugar ring; add all.
                for ring_idx in atom_to_sugar_ring[n_idx]:
                    neighbor_sugar_rings.add(ring_idx)
        if len(neighbor_sugar_rings) >= 2:
            glycosidic_linkage_found = True
            break
    
    if not glycosidic_linkage_found:
        return False, "Candidate sugar rings detected but no glycosidic linkage found connecting them"
    
    return True, f"Found {len(sugar_rings)} candidate sugar ring(s) and at least one glycosidic linkage connecting them"

# Example test cases (the examples provided in the task are very complex).
if __name__ == "__main__":
    # Example SMILES of a disaccharide (alpha-D-Galp-(1->6)-beta-D-Galp) from the provided list.
    test_smiles = "O[C@H]1[C@H](O)O[C@@H]([C@@H]([C@@H]1O)O)CO[C@@H]2[C@@H]([C@H]([C@H]([C@H](O2)CO)O)O)O"
    res, reason = is_oligosaccharide(test_smiles)
    print("Test disaccharide:", res, reason)