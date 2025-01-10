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
    if num_carbons < 25:
        return False, f"Too few carbon atoms ({num_carbons}) to be a triterpenoid"

    # Count number of rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 4:
        return False, f"Too few rings ({num_rings}) to be a triterpenoid"

    # Analyze fused ring systems
    # Get list of rings (atom indices)
    ring_atom_indices = ring_info.AtomRings()
    # Build a list of fused rings
    fused_rings = []
    visited_rings = set()

    for i, ring_i in enumerate(ring_atom_indices):
        if i in visited_rings:
            continue
        fused_ring = set(ring_i)
        visited_rings.add(i)
        # Compare with other rings to see if they are fused (share at least one atom)
        for j, ring_j in enumerate(ring_atom_indices):
            if j <= i or j in visited_rings:
                continue
            if set(ring_i) & set(ring_j):
                fused_ring.update(ring_j)
                visited_rings.add(j)
        fused_rings.append(fused_ring)

    # Find the largest fused ring system
    max_fused_ring_size = max(len(fused_ring) for fused_ring in fused_rings) if fused_rings else 0
    if max_fused_ring_size < 15:
        return False, f"Largest fused ring system contains too few atoms ({max_fused_ring_size}) to be a triterpenoid"

    return True, "Molecule meets the criteria of a triterpenoid"