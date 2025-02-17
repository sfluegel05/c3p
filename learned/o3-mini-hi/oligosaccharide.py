"""
Classifies: CHEBI:50699 oligosaccharide
"""
#!/usr/bin/env python
"""
Classifies: Oligosaccharide 
An oligosaccharide is here defined as a compound formed from at least two monosaccharide units,
each of which is represented by a 5- or 6-membered ring (with one ring oxygen and several hydroxyl substituents)
that are connected by at least one glycosidic (C–O–C) linkage (typically via an exocyclic bridging oxygen).
This revised heuristic:
  1. Adds explicit hydrogens to allow detection of -OH substituents.
  2. Identifies candidate sugar rings:
       • rings of size 5 (furanose) or 6 (pyranose) with exactly one ring oxygen,
         and at least two ring carbons that have an exocyclic -OH (i.e. an oxygen bonded to a hydrogen that is not in the ring).
  3. Merges overlapping rings into distinct sugar units.
  4. Searches for a glycosidic linkage defined as a non‐ring oxygen atom bonded to two carbon atoms, 
     each belonging to a different sugar unit.
If at least two distinct sugar units are found and at least one glycosidic bridge is detected,
the molecule is classified as an oligosaccharide.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_oligosaccharide(smiles: str):
    """
    Determines whether a molecule is an oligosaccharide based on its SMILES string.
    
    Args:
       smiles (str): The SMILES string of the molecule.
    Returns:
       bool: True if the molecule is classified as an oligosaccharide, False otherwise.
       str: Explanation of the classification decision.
    """
    # Parse the molecule and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    candidate_rings = []  # List of sets of atom indices for rings that look like sugar rings.
    # Pre-compile a SMARTS query for a hydroxyl group.
    oh_smarts = Chem.MolFromSmarts("[OX2H]")
    
    # Step 1: Identify candidate sugar rings.
    for ring in atom_rings:
        # Consider only rings of size 5 or 6.
        if len(ring) not in (5, 6):
            continue
        # Count oxygen atoms in the ring.
        ring_oxygen_count = 0
        skip_ring = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Skip rings containing any aromatic atoms.
            if atom.GetIsAromatic():
                skip_ring = True
                break
            if atom.GetAtomicNum() == 8:
                ring_oxygen_count += 1
        if skip_ring:
            continue
        # Expect exactly one ring oxygen.
        if ring_oxygen_count != 1:
            continue
            
        # Next, check for at least two ring carbons bearing an exocyclic hydroxyl.
        oh_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue  # Only consider carbons
            # Look at neighbors that are not in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Check if neighbor matches an -OH pattern.
                if nbr.HasQuery(False): 
                    # In some cases query info is set; here we simply check atomic numbers.
                    continue
                # We require neighbor to be oxygen and have at least one hydrogen.
                if nbr.GetAtomicNum() == 8:
                    # Count attached hydrogens:
                    nHs = sum(1 for n in nbr.GetNeighbors() if n.GetAtomicNum() == 1)
                    if nHs >= 1:
                        oh_count += 1
                        break  # Only count one -OH per carbon
        if oh_count < 2:
            continue
            
        candidate_rings.append(set(ring))
    
    if len(candidate_rings) < 2:
        return False, f"Found {len(candidate_rings)} candidate sugar ring(s); need at least 2 for an oligosaccharide"
    
    # Step 2: Merge overlapping candidate rings into distinct sugar units.
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
    
    # Map each atom index (that belongs to a candidate ring) to its sugar unit id.
    atom_to_unit = {}
    for unit_id, unit_atoms in enumerate(sugar_unit_list):
        for idx in unit_atoms:
            atom_to_unit[idx] = unit_id
            
    # Step 3: Identify glycosidic linkages.
    # Look for non-ring oxygen atoms (exocyclic) that are bonded to carbons belonging to different sugar units.
    glycosidic_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        # Only consider oxygen atoms that are NOT part of any candidate sugar ring,
        # or if they have additional bonds beyond the ring.
        in_ring = False
        for idx in atom.GetIdx(),:
            if idx in atom_to_unit:
                in_ring = True
                break
        # We still accept an oxygen if it is exocyclic (has a neighbor not in the ring)
        if atom.IsInRing():
            # To be conservative, skip oxygens that are only part of a ring.
            continue
            
        connected_units = set()
        # Look at neighbor atoms; if carbon and if that carbon belongs to a sugar unit, record that.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                idx = nbr.GetIdx()
                if idx in atom_to_unit:
                    connected_units.add(atom_to_unit[idx])
        # If the oxygen connects carbons from at least two different sugar units, we count it as a glycosidic linkage.
        if len(connected_units) >= 2:
            glycosidic_found = True
            break
            
    if not glycosidic_found:
        return False, "Candidate sugar units detected but no glycosidic (bridging) bond found connecting distinct units"
    
    return True, f"Found {len(sugar_unit_list)} distinct sugar unit(s) and at least one glycosidic linkage connecting them"

# Example test block.
if __name__ == "__main__":
    # Example: testing on a known disaccharide (alpha-D-Galp-(1->6)-beta-D-Galp).
    test_smiles = "O[C@H]1[C@H](O)O[C@@H]([C@@H]([C@@H]1O)O)CO[C@@H]2[C@@H]([C@H]([C@H]([C@H](O2)CO)O)O)O"
    res, reason = is_oligosaccharide(test_smiles)
    print("Test disaccharide:", res, reason)