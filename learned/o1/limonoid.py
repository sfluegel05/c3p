"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    A limonoid is a triterpenoid that is highly oxygenated and has a prototypical structure
    either containing or derived from a precursor with a 4,4,8-trimethyl-17-furanylsteroid skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Step 1: Check for the steroid nucleus (4 fused rings with sizes 6,6,6,5)
    # Get all rings
    sssr = [list(x) for x in Chem.GetSymmSSSR(mol)]
    if len(sssr) == 0:
        return False, "Molecule has no rings"

    # Map atoms to rings
    from collections import defaultdict
    atom_ring_mapping = defaultdict(list)
    for ring_idx, ring_atoms in enumerate(sssr):
        for atom_idx in ring_atoms:
            atom_ring_mapping[atom_idx].append(ring_idx)

    # Build ring adjacency matrix
    num_rings = len(sssr)
    ring_adj = [[0]*num_rings for _ in range(num_rings)]
    for atom_idx, rings in atom_ring_mapping.items():
        if len(rings) > 1:
            for i in range(len(rings)):
                for j in range(i+1, len(rings)):
                    ring_adj[rings[i]][rings[j]] = 1
                    ring_adj[rings[j]][rings[i]] = 1

    # Find connected ring systems
    def dfs(ring_idx, visited, current_system):
        visited.add(ring_idx)
        current_system.append(ring_idx)
        for neighbor_idx, connected in enumerate(ring_adj[ring_idx]):
            if connected and neighbor_idx not in visited:
                dfs(neighbor_idx, visited, current_system)

    ring_systems = []
    visited_rings = set()
    for ring_idx in range(num_rings):
        if ring_idx not in visited_rings:
            current_system = []
            dfs(ring_idx, visited_rings, current_system)
            ring_systems.append(current_system)

    # Check for steroid nucleus
    steroid_system_found = False
    for system in ring_systems:
        if len(system) == 4:
            # Get ring sizes
            ring_sizes = [len(sssr[ring_idx]) for ring_idx in system]
            ring_sizes.sort()
            if ring_sizes == [5,6,6,6]:
                steroid_system_found = True
                steroid_system = system
                break
    if not steroid_system_found:
        return False, "No steroid nucleus (4 fused rings with sizes 6,6,6,5) found"

    # Step 2: Check for furan ring
    furan_pattern = Chem.MolFromSmarts('c1ccoc1')  # furan ring
    furan_matches = mol.GetSubstructMatches(furan_pattern)
    if not furan_matches:
        return False, "No furan ring found"

    # Step 3: Check if furan ring is connected to steroid nucleus
    # Get atoms in steroid system
    steroid_atom_idxs = set()
    for ring_idx in steroid_system:
        steroid_atom_idxs.update(sssr[ring_idx])

    furan_connected = False
    for furan_match in furan_matches:
        # Check if any atom in furan ring is connected to any atom in steroid system
        for idx in furan_match:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in steroid_atom_idxs:
                    furan_connected = True
                    break
            if furan_connected:
                break
        if furan_connected:
            break
    if not furan_connected:
        return False, "Furan ring not connected to steroid nucleus"

    # Step 4: Check if molecule is highly oxygenated
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, f"Not highly oxygenated, contains only {o_count} oxygen atoms"

    return True, "Molecule is a limonoid (contains steroid nucleus with attached furan ring and is highly oxygenated)"