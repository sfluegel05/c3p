"""
Classifies: CHEBI:36615 triterpenoid
"""
from rdkit import Chem

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    A triterpenoid is derived from a triterpene, having a C30 skeleton which may be rearranged
    or missing some methyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a triterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 24:
        return False, f"Too few carbon atoms ({num_carbons}) to be a triterpenoid"

    # Count number of rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 3:
        return False, f"Too few rings ({num_rings}) to be a triterpenoid"

    # Analyze fused ring systems
    # Get list of rings (as sets of atom indices)
    ring_atom_indices = ring_info.AtomRings()
    rings = [set(ring) for ring in ring_atom_indices]

    # Build adjacency list for rings
    adjacency = {i: set() for i in range(len(rings))}
    for i in range(len(rings)):
        for j in range(i+1, len(rings)):
            if rings[i] & rings[j]:
                adjacency[i].add(j)
                adjacency[j].add(i)

    # Find fused ring systems using BFS
    fused_ring_systems = []
    visited = set()
    for i in range(len(rings)):
        if i not in visited:
            queue = [i]
            visited.add(i)
            fused_system = set(rings[i])
            while queue:
                current = queue.pop(0)
                for neighbor in adjacency[current]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)
                        fused_system.update(rings[neighbor])
            fused_ring_systems.append(fused_system)

    # Find the largest fused ring system
    max_fused_ring_size = max(len(system) for system in fused_ring_systems) if fused_ring_systems else 0
    if max_fused_ring_size < 10:
        return False, f"Largest fused ring system contains too few atoms ({max_fused_ring_size}) to be a triterpenoid"

    return True, "Molecule meets the criteria of a triterpenoid"