"""
Classifies: CHEBI:50699 oligosaccharide
"""
#!/usr/bin/env python
"""
Classifies: Oligosaccharide

An oligosaccharide is defined here as a compound in which at least two monosaccharide units (candidate sugar rings)
are joined by at least one glycosidic (C–O–C) bridging linkage. 

Heuristic improvements over our previous attempt:
  1. We parse the molecule from its SMILES string and add explicit hydrogens.
  2. We scan each ring (using RDKit ring info) and require that candidate rings have size 5 or 6,
     exactly one ring oxygen, and at least a minimal number of exocyclic hydroxyl groups (i.e. an oxygen not in the ring
     that bears at least one hydrogen) attached to ring carbons.
  3. We merge any overlapping candidate rings into one sugar “unit” via a union-find procedure.
  4. We attempt to detect a glycosidic linkage by scanning for oxygen atoms that, via a C–O–C bond,
     connect carbons that are assigned to two distinct candidate sugar units.
If at least two distinct sugar units are found and a glycosidic linkage is detected, the molecule is classified as an oligosaccharide.
Otherwise it is not.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_oligosaccharide(smiles: str):
    """
    Determines whether a molecule is an oligosaccharide based on its SMILES string.
    
    Args:
       smiles (str): SMILES string of the molecule.
    Returns:
       bool: True if classified as an oligosaccharide, else False.
       str: Explanation (reason) for the classification decision.
    """
    # Parse molecule and add hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Identify candidate sugar rings (only rings of size 5 or 6)
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    candidate_rings = []  # Each candidate is a set of atom indices.
    
    for ring in atom_rings:
        if len(ring) not in (5, 6):
            continue
        # Count oxygen atoms in the ring.
        ring_oxy = 0
        skip_ring = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # If any ring atom is aromatic, skip ring (we assume sugar rings are aliphatic)
            if atom.GetIsAromatic():
                skip_ring = True
                break
            if atom.GetAtomicNum() == 8:
                ring_oxy += 1
        if skip_ring:
            continue
        # Require exactly one oxygen in the ring (typical for furanose and pyranose).
        if ring_oxy != 1:
            continue
        
        # Count the number of ring carbons with an exocyclic hydroxyl.
        oh_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # We only consider carbons.
            if atom.GetAtomicNum() != 6:
                continue
            # Check neighbors not in ring for an oxygen that bears at least one hydrogen.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # Count explicit hydrogens on this neighbor.
                    nHs = sum(1 for n in nbr.GetNeighbors() if n.GetAtomicNum() == 1)
                    if nHs >= 1:
                        oh_count += 1
                        break   # Only count once per carbon
        # Here we require at least one or two exocyclic hydroxyls.
        if oh_count < 1:
            continue
        candidate_rings.append(set(ring))
    
    if len(candidate_rings) < 2:
        return False, f"Found {len(candidate_rings)} candidate sugar ring(s); need at least 2 for an oligosaccharide."
    
    # Merge overlapping candidate rings to form distinct sugar units using union-find.
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
                union(i,j)
    sugar_units = {}
    for i, ring in enumerate(candidate_rings):
        root = find(i)
        sugar_units.setdefault(root, set()).update(ring)
    sugar_unit_list = list(sugar_units.values())
    
    if len(sugar_unit_list) < 2:
        return False, "Candidate sugar rings merged into one unit; need at least two distinct monosaccharide units."
    
    # Map each atom (that belongs to any candidate ring) to its sugar unit id.
    atom_to_unit = {}
    for unit_id, unit_atoms in enumerate(sugar_unit_list):
        for idx in unit_atoms:
            atom_to_unit[idx] = unit_id
    
    # Identify glycosidic linkages.
    # We scan for bonds through an oxygen atom that is outside one (or not strictly confined to) candidate sugar rings
    # and that connects two carbons that belong to different sugar units.
    glycosidic_found = False
    for bond in mol.GetBonds():
        # We first check if the bond is a single bond connecting an oxygen and a carbon.
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Order the atoms: let one be oxygen and one be carbon.
        if (a1.GetAtomicNum() == 8 and a2.GetAtomicNum() == 6) or (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 8):
            # Identify the oxygen in the bond.
            oxy_atom = a1 if a1.GetAtomicNum() == 8 else a2
            # We want the oxygen to serve as a bridging atom – it should be connected to at least 2 carbons.
            carbons = [nbr for nbr in oxy_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
            if len(carbons) < 2:
                continue
            # If the oxygen atom is in a ring and all of its bonded atoms are in candidate sugar rings,
            # then likely it is the ring oxygen rather than a bridging (exocyclic) oxygen.
            if oxy_atom.IsInRing():
                extra = any(nbr.GetIdx() not in atom_to_unit for nbr in oxy_atom.GetNeighbors())
                if not extra:
                    continue
            # Determine sugar unit membership for each connected carbon.
            connected_units = set()
            for c in carbons:
                if c.GetIdx() in atom_to_unit:
                    connected_units.add(atom_to_unit[c.GetIdx()])
            # If the oxygen connects carbons from at least 2 different sugar units, mark it as glycosidic.
            if len(connected_units) >= 2:
                glycosidic_found = True
                break

    if not glycosidic_found:
        return False, "Candidate sugar units detected but no glycosidic (bridging) C–O–C bond connecting distinct units was found."
    
    return True, f"Found {len(sugar_unit_list)} distinct sugar unit(s) and at least one glycosidic linkage connecting them."

# Example test
if __name__ == "__main__":
    # Test example: disaccharide alpha-D-Galp-(1->6)-beta-D-Galp.
    test_smiles = "O[C@H]1[C@H](O)O[C@@H]([C@@H]([C@@H]1O)O)CO[C@@H]2[C@@H]([C@H]([C@H]([C@H](O2)CO)O)O)O"
    res, reason = is_oligosaccharide(test_smiles)
    print("Test disaccharide:", res, reason)