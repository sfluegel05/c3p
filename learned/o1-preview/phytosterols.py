"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: phytosterols
"""
from rdkit import Chem

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are sterols similar to cholesterol which occur in plants and 
    vary only in carbon side chains and/or presence or absence of a double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    # Build ring adjacency list to identify fused rings
    ring_adj = [set() for _ in rings]
    for i, ring1 in enumerate(rings):
        for j, ring2 in enumerate(rings):
            if i >= j:
                continue
            if set(ring1) & set(ring2):
                # Rings i and j are fused (share at least one atom)
                ring_adj[i].add(j)
                ring_adj[j].add(i)

    # Find fused ring systems using DFS
    visited = set()
    fused_ring_systems = []
    for i in range(len(rings)):
        if i in visited:
            continue
        stack = [i]
        fused_rings = set()
        while stack:
            ring_idx = stack.pop()
            if ring_idx in visited:
                continue
            visited.add(ring_idx)
            fused_rings.add(ring_idx)
            # Add neighboring rings that are fused
            stack.extend(ring_adj[ring_idx] - visited)
        fused_ring_systems.append(fused_rings)

    # Look for steroid backbone: a fused ring system with 4 rings,
    # consisting of 3 six-membered rings and 1 five-membered ring
    has_steroid_backbone = False
    steroid_ring_system = None
    for ring_system in fused_ring_systems:
        if len(ring_system) != 4:
            continue
        ring_sizes = [len(rings[ring_idx]) for ring_idx in ring_system]
        num_six = ring_sizes.count(6)
        num_five = ring_sizes.count(5)
        if num_six == 3 and num_five == 1:
            has_steroid_backbone = True
            steroid_ring_system = ring_system
            break

    if not has_steroid_backbone:
        return False, "No steroid backbone found"

    # Identify atoms in the steroid backbone
    backbone_atoms = set()
    for ring_idx in steroid_ring_system:
        backbone_atoms.update(rings[ring_idx])

    # Identify side chains attached to the backbone
    side_chain_atoms = set()
    atoms_to_visit = set()
    for atom_idx in backbone_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            nbr_idx = neighbor.GetIdx()
            if nbr_idx not in backbone_atoms and nbr_idx not in side_chain_atoms:
                side_chain_atoms.add(nbr_idx)
                atoms_to_visit.add(nbr_idx)

    # Traverse side chain atoms to get full side chain
    visited_atoms = set()
    while atoms_to_visit:
        current_idx = atoms_to_visit.pop()
        if current_idx in visited_atoms:
            continue
        visited_atoms.add(current_idx)
        current_atom = mol.GetAtomWithIdx(current_idx)
        for neighbor in current_atom.GetNeighbors():
            nbr_idx = neighbor.GetIdx()
            if nbr_idx not in backbone_atoms and nbr_idx not in visited_atoms:
                side_chain_atoms.add(nbr_idx)
                atoms_to_visit.add(nbr_idx)

    # Count the number of carbon atoms in side chains
    num_side_chain_carbons = sum(1 for idx in side_chain_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)

    # Cholesterol has 8 carbons in the side chain; phytosterols typically have more
    if num_side_chain_carbons <= 8:
        return False, f"Side chain too short ({num_side_chain_carbons} carbons), possibly cholesterol"

    return True, "Molecule is a phytosterol with steroid backbone and appropriate side chain length"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'phytosterols',
        'definition': 'Sterols similar to cholesterol which occur in plants and vary only in carbon side chains and/or presence or absence of a double bond.',
        'parents': []
    },
    'config': {
        # Configuration details can be added here if needed
    },
    'message': None,
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    # Additional metadata fields can be added here
}