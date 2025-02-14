"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
"""
Classifies: 3-oxo-Delta(4) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    A 3-oxo-Delta(4) steroid is defined as a steroid with a ketone group at position 3
    and a double bond between carbons 4 and 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(4) steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()

    # Find fused ring systems
    # Initialize a list to hold ring systems (each as a set of ring indices)
    ring_systems = []
    visited = set()
    for i, ring1 in enumerate(atom_rings):
        if i in visited:
            continue
        fused_system = set([i])
        to_visit = [i]
        while to_visit:
            idx = to_visit.pop()
            visited.add(idx)
            ring_atoms = set(atom_rings[idx])
            # Check for fusion with other rings
            for j, ring2 in enumerate(atom_rings):
                if j in visited or j in fused_system:
                    continue
                if ring_atoms & set(atom_rings[j]):
                    fused_system.add(j)
                    to_visit.append(j)
        ring_systems.append(fused_system)

    # Search for a fused ring system with 4 rings of sizes 6,6,6,5
    steroid_system = None
    for system in ring_systems:
        ring_sizes = [len(atom_rings[idx]) for idx in system]
        if sorted(ring_sizes) == [5,6,6,6]:
            steroid_system = system
            break

    if not steroid_system:
        return False, "No steroid nucleus (6-6-6-5 fused ring system) found"

    # Create a substructure of the steroid nucleus
    steroid_atoms = set()
    for idx in steroid_system:
        steroid_atoms.update(atom_rings[idx])
    steroid_submol = Chem.PathToSubmol(mol, list(steroid_atoms))

    # Define SMARTS pattern for ketone at position 3
    ketone_pattern = Chem.MolFromSmarts('O=C[C;R]')
    if not steroid_submol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group at position 3 found"

    # Define SMARTS pattern for double bond between carbons 4 and 5
    double_bond_pattern = Chem.MolFromSmarts('C=C[C;R]')
    if not steroid_submol.HasSubstructMatch(double_bond_pattern):
        return False, "No double bond between carbons 4 and 5 found"

    return True, "Contains steroid nucleus with ketone at position 3 and Δ⁴ double bond"