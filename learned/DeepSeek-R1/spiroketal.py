"""
Classifies: CHEBI:72600 spiroketal
"""
"""
Classifies: spiroketal (CHEBI:71629)
A cyclic ketal in which the ketal carbon is the only common atom of two rings.
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal has a carbon atom shared between two rings (spiro atom),
    each connected via oxygen atoms in a ketal configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a spiroketal, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Get ring information
    ri = mol.GetRingInfo()
    ri.Initialize(mol)  # Ensure ring info is computed
    ring_atoms = ri.AtomRings()

    # Find spiro candidates: atoms in exactly two rings that share no other atoms
    spiro_candidates = []
    for atom in mol.GetAtoms():
        atom_idx = atom.GetIdx()
        rings = ri.AtomRings(atom_idx)
        if len(rings) != 2:
            continue
        # Get the two rings' atom indices
        ring1 = set(ring_atoms[rings[0]])
        ring2 = set(ring_atoms[rings[1]])
        intersection = ring1 & ring2
        if len(intersection) == 1 and atom_idx in intersection:
            spiro_candidates.append(atom)

    # Check each candidate for spiroketal conditions
    for candidate in spiro_candidates:
        # Check for exactly two oxygen neighbors
        neighbors = candidate.GetNeighbors()
        oxygen_neighbors = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(oxygen_neighbors) != 2:
            continue

        # Get the two rings the candidate is part of
        rings = ri.AtomRings(candidate.GetIdx())
        ring1_idx, ring2_idx = rings
        ring1_atoms = set(ring_atoms[ring1_idx])
        ring2_atoms = set(ring_atoms[ring2_idx])

        # Check each oxygen is in a different ring
        o1, o2 = oxygen_neighbors
        o1_in_ring1 = o1.GetIdx() in ring1_atoms
        o1_in_ring2 = o1.GetIdx() in ring2_atoms
        o2_in_ring1 = o2.GetIdx() in ring1_atoms
        o2_in_ring2 = o2.GetIdx() in ring2_atoms

        # Ensure each oxygen is in a different ring
        if (o1_in_ring1 and o2_in_ring2) or (o1_in_ring2 and o2_in_ring1):
            return True, "Spiro atom with two oxygen atoms in separate rings"

    return False, "No spiroketal structure detected"