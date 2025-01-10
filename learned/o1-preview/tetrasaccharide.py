"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide is an oligosaccharide comprising four monomeric monosaccharide units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify rings of size 5 or 6 with exactly one oxygen (monosaccharide rings)
    ring_info = mol.GetRingInfo()
    ring_atoms = ring_info.AtomRings()
    mono_rings = []
    for ring in ring_atoms:
        if len(ring) == 5 or len(ring) == 6:
            o_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if o_count == 1:
                mono_rings.append(ring)

    total_monosaccharides = len(mono_rings)

    if total_monosaccharides != 4:
        return False, f"Found {total_monosaccharides} monosaccharide units, expected 4"

    # Map atoms to monosaccharide units
    atom_to_ring = {}
    for i, ring in enumerate(mono_rings):
        for idx in ring:
            atom_to_ring[idx] = i

    # Identify glycosidic linkages
    glycosidic_bonds = []
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        idx1 = a1.GetIdx()
        idx2 = a2.GetIdx()
        ring1 = atom_to_ring.get(idx1)
        ring2 = atom_to_ring.get(idx2)
        # Bond between atoms in different rings
        if ring1 is not None and ring2 is not None and ring1 != ring2:
            # Check if bond is via an oxygen atom (glycosidic oxygen)
            if a1.GetAtomicNum() == 8 or a2.GetAtomicNum() == 8:
                glycosidic_bonds.append((ring1, ring2))

    unique_linkages = set()
    for r1, r2 in glycosidic_bonds:
        linkage = tuple(sorted((r1, r2)))
        unique_linkages.add(linkage)

    if len(unique_linkages) != 3:
        return False, f"Found {len(unique_linkages)} glycosidic linkages, expected 3"

    # Check if the monosaccharide units are connected via glycosidic linkages
    from collections import deque

    # Build connectivity graph
    graph = {i: set() for i in range(len(mono_rings))}
    for r1, r2 in unique_linkages:
        graph[r1].add(r2)
        graph[r2].add(r1)

    # Perform BFS to check if all units are connected
    visited = set()
    queue = deque()
    queue.append(0)
    visited.add(0)
    while queue:
        current = queue.popleft()
        for neighbor in graph[current]:
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append(neighbor)

    if len(visited) != 4:
        return False, "Monosaccharide units are not properly connected"

    return True, "Molecule contains 4 monosaccharide units connected via glycosidic linkages"