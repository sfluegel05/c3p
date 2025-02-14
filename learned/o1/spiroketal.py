"""
Classifies: CHEBI:72600 spiroketal
"""
from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal is a cyclic ketal in which the ketal carbon is the only common atom of two rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a spiroketal, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for spiroketal
    # [C;R2;X4] - sp3 carbon in exactly two rings
    # Connected to two oxygens [O;R1], each connected to another carbon [C]
    spiroketal_smarts = '[C;R2;X4]([O;R1;$([O][C])])([O;R1;$([O][C])])'

    pattern = Chem.MolFromSmarts(spiroketal_smarts)
    if pattern is None:
        return False, "Invalid SMARTS pattern"

    # Search for the spiroketal pattern in the molecule
    matches = mol.GetSubstructMatches(pattern)

    if not matches:
        return False, "No spiroketal functionality found"

    # Optionally, we can check if the rings share only the spiro carbon
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    for match in matches:
        spiro_atom_idx = match[0]
        oxygen1_idx = match[1]
        oxygen2_idx = match[2]

        # Get the two rings containing the spiro atom
        rings_with_spiro = [ring for ring in atom_rings if spiro_atom_idx in ring]
        if len(rings_with_spiro) != 2:
            continue  # Spiro atom is not in exactly two rings

        # Check that the rings share only the spiro carbon
        shared_atoms = set(rings_with_spiro[0]) & set(rings_with_spiro[1])
        if len(shared_atoms) != 1:
            continue  # Rings share more than one atom

        # Verify that each oxygen is part of a ring
        oxygen1 = mol.GetAtomWithIdx(oxygen1_idx)
        oxygen2 = mol.GetAtomWithIdx(oxygen2_idx)
        if not oxygen1.IsInRing() or not oxygen2.IsInRing():
            continue  # Oxygen atoms are not in rings

        # If all checks passed, we have found a spiroketal
        return True, f"Molecule is a spiroketal with spiro ketal center at atom index {spiro_atom_idx}"

    return False, "No spiroketal functionality found"