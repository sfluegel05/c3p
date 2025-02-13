"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies: sterol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is any 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol
    (additional carbon atoms may be present in the side chain).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get the ring information
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    # Must have at least 4 rings
    if len(rings) < 4:
        return False, f"Contains only {len(rings)} rings, requires at least 4"

    # Build ring adjacency graph
    ring_count = len(rings)
    ring_adjacency = [set() for _ in range(ring_count)]

    # Rings are connected if they share at least 2 atoms (fused)
    for i in range(ring_count):
        for j in range(i+1, ring_count):
            shared_atoms = set(rings[i]) & set(rings[j])
            if len(shared_atoms) >= 2:
                ring_adjacency[i].add(j)
                ring_adjacency[j].add(i)

    # Find fused ring systems
    visited = set()
    fused_ring_systems = []

    for i in range(ring_count):
        if i not in visited:
            stack = [i]
            component = []
            while stack:
                ring_idx = stack.pop()
                if ring_idx not in visited:
                    visited.add(ring_idx)
                    component.append(ring_idx)
                    stack.extend(ring_adjacency[ring_idx] - visited)
            fused_ring_systems.append(component)

    # Check if there is a fused ring system with at least 4 rings
    steroid_ring_system = None
    for system in fused_ring_systems:
        if len(system) >= 4:
            # Check if ring sizes are correct (three 6-membered and one 5-membered ring)
            ring_sizes = sorted([len(rings[i]) for i in system])
            if ring_sizes[:4] == [5,6,6,6]:
                steroid_ring_system = system
                break

    if steroid_ring_system is None:
        return False, "Does not contain steroid nucleus with fused 5,6,6,6 rings"

    # Define SMARTS pattern for 3-hydroxy group attached to ring A (first ring)
    hydroxy_pattern = Chem.MolFromSmarts('[#6]-1([#8H])-[#6]-[#6]-[#6]-[#6]-1')  # 3-hydroxy group
    if hydroxy_pattern is None:
        return False, "Invalid hydroxy group SMARTS pattern"

    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "Does not have 3-hydroxy group at position 3"

    return True, "Contains steroid nucleus with 3-hydroxy group characteristic of sterols"

__metadata__ = {
    'chemical_class': {
        'name': 'sterol',
        'definition': 'Any 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol (additional carbon atoms may be present in the side chain).'
    }
}