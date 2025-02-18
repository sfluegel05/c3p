"""
Classifies: CHEBI:50699 oligosaccharide
"""
#!/usr/bin/env python
"""
Classifies: Oligosaccharide 
An oligosaccharide is here defined as a compound formed from at least two monosaccharide units,
each represented by a 5- or 6-membered ring (with one ring oxygen and several hydroxyl substituents)
that are connected by at least one glycosidic (C–O–C) linkage (typically via an exocyclic bridging oxygen).

Heuristic:
  1. Parse the molecule from SMILES and add explicit hydrogens.
  2. Identify candidate sugar rings (only rings of size 5 or 6 having exactly one ring oxygen and at least two ring carbons 
     that bear an exocyclic hydroxyl group – i.e. an oxygen atom bonded to a hydrogen that is not part of the ring).
  3. Merge overlapping rings into distinct sugar units.
  4. Identify a glycosidic linkage defined here as a non‐ring oxygen atom connected to two carbons from different sugar units.
If at least two sugar units are detected and a glycosidic linkage is found, the molecule is classified as an oligosaccharide.
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
    
    candidate_rings = []  # List of sets of atom indices that look like sugar rings.
    
    # Identify candidate sugar rings.
    for ring in atom_rings:
        # Consider only rings of size 5 or 6.
        if len(ring) not in (5, 6):
            continue
        # Count oxygen atoms in the ring.
        ring_oxygen_count = 0
        skip_ring = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Skip rings that contain any aromatic atoms.
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
            
        # Check for at least two ring carbons bearing an exocyclic hydroxyl.
        oh_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue  # Only consider carbons
            # Look at neighbors that are not in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Check if neighbor is oxygen and has at least one hydrogen attached.
                if nbr.GetAtomicNum() == 8:
                    nHs = sum(1 for n in nbr.GetNeighbors() if n.GetAtomicNum() == 1)
                    if nHs >= 1:
                        oh_count += 1
                        break  # Only count one -OH per carbon
        if oh_count < 2:
            continue
            
        candidate_rings.append(set(ring))
    
    if len(candidate_rings) < 2:
        return False, f"Found {len(candidate_rings)} candidate sugar ring(s); need at least 2 for an oligosaccharide."
    
    # Merge overlapping candidate rings into distinct sugar units.
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
        return False, "Candidate sugar rings merged into one unit; need at least two distinct monosaccharide units."
    
    # Map each atom index (that belongs to a candidate ring) to its sugar unit id.
    atom_to_unit = {}
    for unit_id, unit_atoms in enumerate(sugar_unit_list):
        for idx in unit_atoms:
            atom_to_unit[idx] = unit_id
            
    # Identify glycosidic linkages.
    # Look for exocyclic oxygen atoms (not part of a candidate sugar ring or only in ring) that connect carbons from different sugar units.
    glycosidic_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        # If an oxygen atom is in a ring and has no extra-bond outside, skip it.
        if atom.IsInRing():
            # Check neighbors: if all neighbors are in a candidate sugar unit
            # then it likely serves only as a ring heteroatom. We require at least one neighbor outside.
            extra_neighbor = False
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in atom_to_unit:
                    extra_neighbor = True
                    break
            if not extra_neighbor:
                continue

        connected_units = set()
        # Look at neighboring atoms; if carbon and if that carbon belongs to a sugar unit, record that unit.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                idx = nbr.GetIdx()
                if idx in atom_to_unit:
                    connected_units.add(atom_to_unit[idx])
        if len(connected_units) >= 2:
            glycosidic_found = True
            break
            
    if not glycosidic_found:
        return False, "Candidate sugar units detected but no glycosidic (bridging) bond was found connecting distinct units."
    
    return True, f"Found {len(sugar_unit_list)} distinct sugar unit(s) and at least one glycosidic linkage connecting them."

# Example test block.
if __name__ == "__main__":
    # Example: testing on a known disaccharide (alpha-D-Galp-(1->6)-beta-D-Galp).
    test_smiles = "O[C@H]1[C@H](O)O[C@@H]([C@@H]([C@@H]1O)O)CO[C@@H]2[C@@H]([C@H]([C@H]([C@H](O2)CO)O)O)O"
    res, reason = is_oligosaccharide(test_smiles)
    print("Test disaccharide:", res, reason)