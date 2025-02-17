"""
Classifies: CHEBI:50699 oligosaccharide
"""
#!/usr/bin/env python
"""
Classifies: Oligosaccharide
An oligosaccharide is defined as a compound in which at least two monosaccharide units 
(typically 5- or 6-membered rings with one ring oxygen and several hydroxyl substituents)
are joined by glycosidic (typically C–O–C) linkages.
This heuristic:
  1. Identifies candidate sugar rings: rings of size 5 (furanose) or 6 (pyranose) that contain exactly one oxygen and no aromatic atoms.
  2. Merges overlapping candidate rings into distinct sugar units.
  3. Searches for a bridging (glycosidic) oxygen that connects two sugar units via C–O–C bonds.
If at least two distinct sugar units are found and at least one glycosidic linkage is detected,
the molecule is tentatively classified as an oligosaccharide.
"""

from rdkit import Chem

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    Steps:
      1. Identify candidate sugar rings (5 or 6 membered rings with exactly 1 oxygen and no aromatic atoms).
      2. Group overlapping candidate rings into distinct sugar units.
      3. Search for a glycosidic linkage defined here as a bridging oxygen atom that connects a carbon 
         in one candidate sugar unit to a carbon in a different candidate sugar unit via a C–O and O–C bond.
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
    candidate_rings = []  # Each candidate ring is stored as a set of atom indices.
    for ring in atom_rings:
        # Only 5- and 6-membered rings are considered.
        if len(ring) not in (5,6):
            continue
        oxygen_count = 0
        skip = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Skip ring if any atom is aromatic.
            if atom.GetIsAromatic():
                skip = True
                break
            if atom.GetAtomicNum() == 8:
                oxygen_count += 1
        if skip:
            continue
        # We expect a candidate sugar ring to have exactly one ring oxygen.
        if oxygen_count != 1:
            continue
        candidate_rings.append(set(ring))
    
    if len(candidate_rings) < 2:
        return False, f"Found {len(candidate_rings)} candidate sugar ring(s); need at least 2 for an oligosaccharide"
    
    # Step 2: Merge overlapping candidate rings into distinct sugar units.
    # We use a simple union-find approach.
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
            
    for i in range(len(candidate_rings)):
        for j in range(i+1, len(candidate_rings)):
            if candidate_rings[i] & candidate_rings[j]:
                union(i, j)
    
    sugar_units = {}
    for i, ring in enumerate(candidate_rings):
        root = find(i)
        sugar_units.setdefault(root, set()).update(ring)
    sugar_unit_list = list(sugar_units.values())
    
    if len(sugar_unit_list) < 2:
        return False, "Candidate sugar rings merged into one unit; need at least two distinct monosaccharide units"
    
    # Build a mapping from atom index to sugar unit id.
    atom_to_unit = {}
    for unit_id, unit_atoms in enumerate(sugar_unit_list):
        for idx in unit_atoms:
            atom_to_unit[idx] = unit_id
    
    # Step 3: Look for a glycosidic linkage.
    # Our approach: look for an oxygen atom (typically bridging) that bonds to two carbons
    # where each carbon belongs to a candidate sugar unit, and those units are distinct.
    glycosidic_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue  # Only consider oxygen atoms as potential bridging atoms.
        # Evaluate neighbors of this oxygen.
        connected_units = set()
        carbon_neighbors = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # Consider only carbon neighbors.
                carbon_neighbors.append(nbr)
                idx = nbr.GetIdx()
                if idx in atom_to_unit:
                    connected_units.add(atom_to_unit[idx])
        # Look for a bridging oxygen: it must be bonded to at least two carbons belonging to different sugar units.
        if len(connected_units) >= 2:
            glycosidic_found = True
            break

    if not glycosidic_found:
        return False, "Candidate sugar units detected but no glycosidic (bridging) bond found connecting distinct units"
    
    return True, f"Found {len(sugar_unit_list)} distinct sugar unit(s) and at least one glycosidic linkage connecting them"

# Example test block.
if __name__ == "__main__":
    # Testing on a known disaccharide: alpha-D-Galp-(1->6)-beta-D-Galp.
    test_smiles = "O[C@H]1[C@H](O)O[C@@H]([C@@H]([C@@H]1O)O)CO[C@@H]2[C@@H]([C@H]([C@H]([C@H](O2)CO)O)O)O"
    res, reason = is_oligosaccharide(test_smiles)
    print("Test disaccharide:", res, reason)