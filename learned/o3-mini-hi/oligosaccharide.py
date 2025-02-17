"""
Classifies: CHEBI:50699 oligosaccharide
"""
#!/usr/bin/env python
"""
Classifies: Oligosaccharide
An oligosaccharide is defined as a compound in which at least two monosaccharide units 
(typically 5- or 6-membered rings having one ring oxygen, with several hydroxyl substituents)
are joined by glycosidic (typically C–O–C) linkages. 
This improved heuristic:
  1. Identifies candidate sugar rings (rings of size 5 or 6 with exactly one oxygen and no aromatic atoms).
  2. Groups overlapping candidate rings into single sugar units (e.g. a fused bicyclic system is one unit).
  3. Scans non‐ring bonds for a C–O–C bond that directly connects atoms that belong to different sugar units.
If at least 2 sugar units are detected and at least one glycosidic linkage is found connecting them,
the molecule is tentatively called an oligosaccharide.
"""

from rdkit import Chem

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    Steps:
      1. Identify candidate sugar rings: rings of size 5 (furanose) or 6 (pyranose) with exactly one oxygen & no aromatic atoms.
      2. Merge overlapping candidate rings into a single sugar unit.
      3. Search for a non‐ring bond (typically a C–O–C) that connects atoms from two distinct sugar units.
    Args:
      smiles (str): SMILES string of the molecule.
    Returns:
      bool: True if the molecule is classified as an oligosaccharide, False otherwise.
      str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # Each ring is a tuple of atom indices.
    
    # Step 1: Identify candidate sugar rings.
    candidate_rings = []   # List of sets of atom indices
    for ring in atom_rings:
        if len(ring) not in (5,6):
            continue
        oxygen_count = 0
        skip = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                oxygen_count += 1
            if atom.GetIsAromatic():
                skip = True
                break
        if skip:
            continue
        # Typical sugar ring has exactly one oxygen.
        if oxygen_count != 1:
            continue
        candidate_rings.append(set(ring))
    
    # If we have fewer than 2 candidate rings, then no oligosaccharide
    if len(candidate_rings) < 2:
        return False, f"Found {len(candidate_rings)} candidate sugar ring(s); need at least 2 for oligosaccharide"
    
    # Step 2: Group candidate rings that overlap (share atoms) 
    # We use a simple union-find structure.
    parent = list(range(len(candidate_rings)))
    
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    
    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[ry] = rx
    
    # Merge candidate rings that share at least one atom.
    for i in range(len(candidate_rings)):
        for j in range(i+1, len(candidate_rings)):
            if candidate_rings[i] & candidate_rings[j]:
                union(i, j)
    
    # Now form sugar_units: dict mapping group root -> union of atom indices from candidate rings
    sugar_units = {}
    for i, ring in enumerate(candidate_rings):
        root = find(i)
        sugar_units.setdefault(root, set()).update(ring)
    sugar_unit_list = list(sugar_units.values())
    
    if len(sugar_unit_list) < 2:
        return False, "Candidate sugar rings merged into one unit; need at least two distinct monosaccharide units"
    
    # To ease lookup, create a mapping: atom idx -> sugar unit id (if atom belongs to more than one unit, pick one)
    atom_to_unit = {}
    for unit_id, unit_atoms in enumerate(sugar_unit_list):
        for idx in unit_atoms:
            atom_to_unit[idx] = unit_id
    
    # Step 3: Search for a glycosidic linkage.
    # A typical glycosidic bond is a non‐ring bond (i.e. not part of any ring) connecting atoms that belong
    # to different sugar units; usually the bond is C-O where the oxygen is bridging.
    glycosidic_found = False
    for bond in mol.GetBonds():
        if bond.IsInRing():
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        unit_ids = set()
        # If an atom is part of a sugar unit then record its unit id.
        if a1.GetIdx() in atom_to_unit:
            unit_ids.add(atom_to_unit[a1.GetIdx()])
        if a2.GetIdx() in atom_to_unit:
            unit_ids.add(atom_to_unit[a2.GetIdx()])
        # Look for bond connecting different sugar units.
        if len(unit_ids) >= 2:
            # Also require that at least one end is an oxygen (as glycosidic bonds are C-O-C)
            if a1.GetAtomicNum() == 8 or a2.GetAtomicNum() == 8:
                glycosidic_found = True
                break
    
    if not glycosidic_found:
        return False, "Candidate sugar units detected but no glycosidic (bridging) bond found connecting distinct units"
    
    return True, f"Found {len(sugar_unit_list)} distinct sugar unit(s) and glycosidic linkage connecting them"

# Example small test. (Note: many of the provided examples are very complex.)
if __name__ == "__main__":
    # Test on a known disaccharide: alpha-D-Galp-(1->6)-beta-D-Galp.
    test_smiles = "O[C@H]1[C@H](O)O[C@@H]([C@@H]([C@@H]1O)O)CO[C@@H]2[C@@H]([C@H]([C@H]([C@H](O2)CO)O)O)O"
    res, reason = is_oligosaccharide(test_smiles)
    print("Test disaccharide:", res, reason)