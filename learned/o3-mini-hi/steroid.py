"""
Classifies: CHEBI:35341 steroid
"""
"""
Classifies: Steroid compounds based on a cyclopenta[a]phenanthrene-derived nucleus.
Definition: Any compound (natural or synthetic) based on the cyclopenta[a]phenanthrene skeleton,
which traditionally has three fused six-membered rings and one five-membered ring, though
modifications (bond scissions, ring expansions/contractions) are allowed.
"""

from rdkit import Chem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    The method first looks for an intact cyclopenta[a]phenanthrene pattern using a simplified
    SMARTS. If not found then it examines the fused ring systems for a cluster of rings that matches
    (heuristically) the expected sizes (one 5-membered and at least two 6-membered rings).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First, try to match an intact steroid nucleus (a simplified cyclopenta[a]phenanthrene).
    # Note: This SMARTS is a heuristic representation and may not cover all steroid variations.
    steroid_smarts = "C1CCC2C3CCC4C(C3)C2CC1C4"
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if steroid_pattern is not None and mol.HasSubstructMatch(steroid_pattern):
        return True, "Contains intact cyclopenta[a]phenanthrene skeleton"

    # If no intact steroid nucleus is found, try to detect a fused ring system characteristic of steroids.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # each ring is a tuple of atom indices
    if not rings:
        return False, "No ring system found"

    # Build a graph where each node is a ring (indexed by its position in 'rings')
    # and an edge connects two rings if they share at least 2 atoms (i.e. fused ring).
    ring_graph = {i: set() for i in range(len(rings))}
    for i in range(len(rings)):
        for j in range(i+1, len(rings)):
            if len(set(rings[i]).intersection(rings[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)

    # Find connected components in the ring graph. Each component is a set of ring indices.
    def dfs(node, visited, component):
        visited.add(node)
        component.add(node)
        for neigh in ring_graph[node]:
            if neigh not in visited:
                dfs(neigh, visited, component)
    
    visited = set()
    components = []
    for node in ring_graph:
        if node not in visited:
            comp = set()
            dfs(node, visited, comp)
            components.append(comp)
    
    # Check each connected component for steroid-like ring sizes.
    # Heuristic: an intact steroid nucleus has 4 rings (three rings of size 6 and one of size 5),
    # but due to modifications we allow a cluster with at least 3 rings containing at least one 5-membered
    # and at least two 6-membered rings.
    for comp in components:
        if len(comp) >= 3:
            size5_count = 0
            size6_count = 0
            for i in comp:
                ring_size = len(rings[i])
                if ring_size == 5:
                    size5_count += 1
                elif ring_size == 6:
                    size6_count += 1
            if size5_count >= 1 and size6_count >= 2:
                return True, "Contains a fused ring system resembling a steroid nucleus (at least one 5-membered and two 6-membered rings fused)"

    return False, "No steroid nucleus pattern detected"