"""
Classifies: CHEBI:33848 polycyclic arene
"""
from rdkit import Chem

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene (polycyclic aromatic hydrocarbon) based on its SMILES string.
    A polycyclic arene is a hydrocarbon consisting of fused aromatic rings.
    The molecule may contain substituents with heteroatoms, but the ring system should consist only of carbon atoms.
    The fused ring system should have at least two rings, and at least two rings should be aromatic.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polycyclic arene, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    num_rings = len(atom_rings)
    if num_rings < 2:
        return False, "Molecule does not have at least two rings"

    # Build ring adjacency graph
    num_atoms = mol.GetNumAtoms()
    atom_ring_indices = [[] for _ in range(num_atoms)]
    for ring_idx, ring in enumerate(atom_rings):
        for atom_idx in ring:
            atom_ring_indices[atom_idx].append(ring_idx)

    ring_graph = {}
    for i in range(num_rings):
        ring_graph[i] = set()

    # Build ring adjacency via shared bonds
    bond_rings = ring_info.BondRings()
    for i in range(num_rings):
        for j in range(i+1, num_rings):
            # Check if rings i and j share at least one bond
            bonds_i = set(bond_rings[i])
            bonds_j = set(bond_rings[j])
            if bonds_i & bonds_j:
                ring_graph[i].add(j)
                ring_graph[j].add(i)

    # Find connected components (fused ring systems)
    visited = [False] * num_rings
    ring_systems = []
    for i in range(num_rings):
        if not visited[i]:
            stack = [i]
            component = []
            while stack:
                node = stack.pop()
                if not visited[node]:
                    visited[node] = True
                    component.append(node)
                    stack.extend(ring_graph[node])
            ring_systems.append(component)

    # For each ring system, check if it meets criteria
    for system in ring_systems:
        if len(system) < 2:
            continue  # Skip ring systems with less than two rings
        # Collect atoms in the ring system
        ring_system_atoms = set()
        for ring_idx in system:
            ring_atoms = atom_rings[ring_idx]
            ring_system_atoms.update(ring_atoms)
        # Check that all atoms in the ring system are carbons
        all_carbons = True
        for atom_idx in ring_system_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:
                all_carbons = False
                break
        if not all_carbons:
            continue  # Skip systems with heteroatoms in rings
        # Check that at least two rings are aromatic
        num_aromatic_rings = 0
        for ring_idx in system:
            ring_atoms = atom_rings[ring_idx]
            # Check if ring is aromatic
            ring_bonds = [mol.GetBondBetweenAtoms(ring_atoms[i], ring_atoms[(i+1)%len(ring_atoms)]) for i in range(len(ring_atoms))]
            if all(bond.GetIsAromatic() for bond in ring_bonds):
                num_aromatic_rings += 1
        if num_aromatic_rings >= 2:
            return True, "Molecule contains fused ring system with at least two aromatic rings made of carbon atoms"

    return False, "No suitable fused aromatic ring system with at least two aromatic rings found"