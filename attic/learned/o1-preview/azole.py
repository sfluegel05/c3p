"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: CHEBI:23125 azole
"""

from rdkit import Chem

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole based on its SMILES string.
    An azole is any monocyclic heteroarene consisting of a five-membered ring containing nitrogen.
    Azoles can also contain one or more other non-carbon atoms, such as nitrogen, sulfur, or oxygen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an azole, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ri = mol.GetRingInfo()

    # Get list of atom indices in rings
    atom_rings = ri.AtomRings()

    # Loop through the rings
    for ring_atoms in atom_rings:
        # Check if ring is 5-membered
        if len(ring_atoms) != 5:
            continue  # Not a five-membered ring
        # Check if ring is aromatic
        if not all(mol.GetAtomWithIdx(atom_idx).GetIsAromatic() for atom_idx in ring_atoms):
            continue  # Ring is not aromatic
        # Check if ring contains nitrogen
        if not any(mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 7 for atom_idx in ring_atoms):
            continue  # Ring does not contain nitrogen
        # Check if ring contains only allowed atoms (C, N, O, S)
        if not all(mol.GetAtomWithIdx(atom_idx).GetAtomicNum() in [6, 7, 8, 16] for atom_idx in ring_atoms):
            continue  # Ring contains other atoms
        # Check if ring is monocyclic (each atom is only in one ring)
        if not all(ri.NumAtomRings(atom_idx) == 1 for atom_idx in ring_atoms):
            continue  # Ring is fused with other rings
        # All criteria met, this is an azole
        return True, "Contains monocyclic five-membered heteroaromatic ring with nitrogen"

    # No azole rings found
    return False, "No azole ring found"