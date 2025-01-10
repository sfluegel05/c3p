"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: hopanoid
"""
from rdkit import Chem

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    A hopanoid is a triterpenoid based on a hopane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string to create an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove hydrogens from the molecule
    mol = Chem.RemoveHs(mol)

    # Get ring information
    ri = mol.GetRingInfo()
    ring_count = ri.NumRings()
    if ring_count < 5:
        return False, f"Contains {ring_count} rings, less than 5 rings required for hopanoid"

    # Get ring atom indices
    rings = ri.AtomRings()
    ring_sizes = [len(r) for r in rings]

    # Count six-membered and five-membered rings
    num_six_membered = ring_sizes.count(6)
    num_five_membered = ring_sizes.count(5)

    if num_six_membered < 4 or num_five_membered < 1:
        return False, f"Requires at least 4 six-membered rings and 1 five-membered ring (found {num_six_membered} six-membered and {num_five_membered} five-membered rings)"

    # Check for fused rings
    # Build an adjacency list of rings
    ring_adjacency = {i: set() for i in range(len(rings))}

    for i in range(len(rings)):
        for j in range(i+1, len(rings)):
            # Check if rings i and j share at least two atoms (fused)
            shared_atoms = set(rings[i]) & set(rings[j])
            if len(shared_atoms) >= 2:
                ring_adjacency[i].add(j)
                ring_adjacency[j].add(i)

    # Find connected components (fused ring systems)
    def dfs(ring_idx, visited):
        visited.add(ring_idx)
        for neighbor in ring_adjacency[ring_idx]:
            if neighbor not in visited:
                dfs(neighbor, visited)

    # Find all fused ring systems
    fused_systems = []
    visited_rings = set()
    for i in range(len(rings)):
        if i not in visited_rings:
            component = set()
            dfs(i, component)
            fused_systems.append(component)
            visited_rings.update(component)

    # Look for a fused ring system with at least 5 rings
    hopanoid_system = None
    for system in fused_systems:
        if len(system) >= 5:
            # Check if the system contains at least 4 six-membered and 1 five-membered rings
            system_ring_sizes = [ring_sizes[i] for i in system]
            six_membered = system_ring_sizes.count(6)
            five_membered = system_ring_sizes.count(5)
            if six_membered >= 4 and five_membered >= 1:
                hopanoid_system = system
                break

    if hopanoid_system is not None:
        return True, "Contains ring system characteristic of hopanoids (at least 5 fused rings with 4 six-membered and 1 five-membered ring)"
    else:
        return False, "Does not have the fused ring system characteristic of hopanoids"