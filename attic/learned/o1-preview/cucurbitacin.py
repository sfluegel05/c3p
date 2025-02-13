"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Classifies: cucurbitacin
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids derived from cucurbitane.
    They have a characteristic fused ring system of three six-membered rings
    and one five-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # List of tuples of atom indices
    num_rings = len(atom_rings)
    if num_rings < 4:
        return False, f"Contains {num_rings} rings, less than 4 required for tetracyclic system"

    # Build a list of ring atoms and ring sizes
    ring_atoms_list = [set(ring) for ring in atom_rings]
    ring_sizes = [len(ring) for ring in atom_rings]

    # Build ring fusion graph
    # Nodes are ring indices, edges exist if rings share two or more atoms
    ring_graph = {}
    for i in range(num_rings):
        ring_graph[i] = set()
        for j in range(num_rings):
            if i != j:
                # Check if rings i and j are fused
                shared_atoms = ring_atoms_list[i] & ring_atoms_list[j]
                if len(shared_atoms) >= 2:
                    ring_graph[i].add(j)

    # Find connected components (fused ring systems)
    visited = set()
    fused_ring_systems = []

    for i in range(num_rings):
        if i not in visited:
            # Perform DFS to find all rings in this fused system
            stack = [i]
            component = set()
            while stack:
                ring_idx = stack.pop()
                if ring_idx not in visited:
                    visited.add(ring_idx)
                    component.add(ring_idx)
                    stack.extend(ring_graph[ring_idx] - visited)
            fused_ring_systems.append(component)

    # Look for a fused ring system with four rings: 3 six-membered and 1 five-membered
    tetracyclic_found = False
    for system in fused_ring_systems:
        if len(system) == 4:
            sizes = [ring_sizes[idx] for idx in system]
            num_6 = sizes.count(6)
            num_5 = sizes.count(5)
            if num_6 == 3 and num_5 == 1:
                tetracyclic_found = True
                break

    if not tetracyclic_found:
        return False, "Does not contain the characteristic tetracyclic ring system of cucurbitacins"

    # Check for triterpenoid skeleton (30 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30:
        return False, f"Carbon count {c_count} is less than 30, not a triterpenoid"

    # Check for characteristic functional groups
    # Cucurbitacins often have multiple hydroxyl and ketone groups
    num_hydroxyl = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[OX2H]')))
    num_ketone = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[$([CX3]=O)]')))
    if num_hydroxyl + num_ketone < 3:
        return False, f"Contains {num_hydroxyl} hydroxyl and {num_ketone} ketone groups, less than 3 total"

    # Passed all checks
    return True, "Contains characteristic tetracyclic ring system and functional groups of cucurbitacins"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'cucurbitacin',
        'definition': 'Any one of a class of tetracyclic triterpenoids, formally derived from the triterpene hydrocarbon cucurbitane, developed by some plants (especially those of the family Cucurbitaceaeas) as a defence mechanism against herbivores.',
        'parents': []
    }
}