"""
Classifies: CHEBI:38077 polypyrrole
"""
"""
Classifies: polypyrrole
"""
from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is a compound composed of two or more pyrrole units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypyrrole, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    pyrrole_count = 0

    # Iterate over rings
    for ring in atom_rings:
        # Check if ring has 5 atoms
        if len(ring) != 5:
            continue

        # Check if ring is aromatic
        is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        if not is_aromatic:
            continue

        # Count number of nitrogen atoms in the ring
        n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)

        if n_count == 1:
            # Found a pyrrole unit
            pyrrole_count += 1

    if pyrrole_count >= 2:
        return True, f"Contains {pyrrole_count} pyrrole units"
    else:
        return False, f"Contains only {pyrrole_count} pyrrole unit(s), less than 2 required for polypyrrole"