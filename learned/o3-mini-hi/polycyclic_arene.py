"""
Classifies: CHEBI:33848 polycyclic arene
"""
#!/usr/bin/env python
"""
Classifies: Polycyclic Arene, defined as a polycyclic aromatic hydrocarbon.
A polycyclic arene is characterized by having two or more fused aromatic rings.
Fused rings share at least one bond (thus at least two atoms).
"""

from rdkit import Chem

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene (polycyclic aromatic hydrocarbon).
    
    The criteria used:
      - SMILES string must be valid.
      - The molecule must contain ring structures.
      - At least two rings must be fully aromatic.
      - At least one set of these aromatic rings is fused (sharing at least two atoms, i.e. a bond).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a polycyclic arene, False otherwise.
        str: Explanation of the classification.
    """
    # Parse SMILES into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # tuple of tuples of atom indices making each ring
    if not atom_rings:
        return False, "No rings detected in the molecule"
    
    # Identify aromatic rings (only consider rings in which every atom is aromatic)
    aromatic_rings = []
    for ring in atom_rings:
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_rings.append(set(ring))
    
    if len(aromatic_rings) < 2:
        return False, "Fewer than two aromatic rings detected"
    
    # Build a graph where each node corresponds to an aromatic ring.
    # Two rings are considered "fused" if they share at least 2 atoms (i.e., a common bond).
    n = len(aromatic_rings)
    ring_graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            # Intersection giving common atoms between rings
            if len(aromatic_rings[i].intersection(aromatic_rings[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Perform depth-first search (DFS) on the ring graph to locate connected components
    visited = set()
    connected_components = []
    for i in range(n):
        if i not in visited:
            stack = [i]
            comp = set()
            while stack:
                node = stack.pop()
                if node not in comp:
                    comp.add(node)
                    stack.extend(ring_graph[node] - comp)
            visited.update(comp)
            connected_components.append(comp)
    
    # If any connected component has 2 or more fused aromatic rings, we have a polycyclic arene
    for comp in connected_components:
        if len(comp) >= 2:
            return True, f"Found a fused aromatic system with {len(comp)} rings"
    
    return False, "No fused aromatic rings (polycyclic system) detected"