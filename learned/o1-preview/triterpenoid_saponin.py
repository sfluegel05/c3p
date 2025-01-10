"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
from rdkit import Chem

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    A triterpenoid saponin is a terpene glycoside in which the terpene moiety is a triterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    bond_rings = ring_info.BondRings()
    if len(bond_rings) == 0:
        return False, "No rings found in the molecule"

    # Build ring adjacency graph
    ring_bond_sets = [set(ring) for ring in bond_rings]
    num_rings = len(ring_bond_sets)
    adj_matrix = [[] for _ in range(num_rings)]

    for i in range(num_rings):
        for j in range(i+1, num_rings):
            # If rings share bonds, they are fused
            if ring_bond_sets[i] & ring_bond_sets[j]:
                # Add edge between rings i and j
                adj_matrix[i].append(j)
                adj_matrix[j].append(i)

    # Find connected components in the ring adjacency graph
    visited = [False] * num_rings
    fused_ring_systems = []

    def dfs(i, component):
        visited[i] = True
        component.append(i)
        for neighbor in adj_matrix[i]:
            if not visited[neighbor]:
                dfs(neighbor, component)

    for i in range(num_rings):
        if not visited[i]:
            component = []
            dfs(i, component)
            fused_ring_systems.append(component)

    # Check for fused ring system with at least five rings
    has_pentacyclic_core = False
    pentacyclic_core_rings = []
    for system in fused_ring_systems:
        if len(system) >= 5:
            has_pentacyclic_core = True
            pentacyclic_core_rings = system
            break

    if not has_pentacyclic_core:
        return False, "No pentacyclic fused ring system found"

    # Look for attached sugar moieties
    # Define a SMARTS pattern for pyranose (six-membered sugar ring)
    sugar_smarts = '[#6&R1]-1-[#6&R1]-[#6&R1]-[#6&R1]-[#6&R1]-[#8&R1]-1'
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    if sugar_pattern is None:
        return False, "Failed to parse sugar SMARTS pattern"

    # Find sugar rings
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) == 0:
        return False, "No sugar moiety found attached"

    # Get atoms in pentacyclic core
    core_atoms = set()
    for ring_idx in pentacyclic_core_rings:
        ring_bonds = bond_rings[ring_idx]
        for bond_idx in ring_bonds:
            bond = mol.GetBondWithIdx(bond_idx)
            core_atoms.add(bond.GetBeginAtomIdx())
            core_atoms.add(bond.GetEndAtomIdx())

    # Check for bonds between core atoms and sugar atoms
    found_linkage = False
    for match in sugar_matches:
        sugar_atoms = set(match)
        # Check if any sugar atom is directly connected to a core atom via an oxygen
        for sugar_atom_idx in sugar_atoms:
            sugar_atom = mol.GetAtomWithIdx(sugar_atom_idx)
            if sugar_atom.GetAtomicNum() != 8:  # We're looking for oxygen atoms
                continue
            for neighbor in sugar_atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx in core_atoms:
                    # Found an oxygen linking sugar to core
                    found_linkage = True
                    break
            if found_linkage:
                break
        if found_linkage:
            break

    if not found_linkage:
        return False, "No glycosidic linkage between core and sugar found"

    return True, "Contains pentacyclic triterpenoid core with attached sugar moieties via glycosidic bonds"