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
    A polypyrrole is defined as a compound composed of two or more pyrrole units.
    Pyrrole units are five-membered rings containing exactly one nitrogen atom.

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

    # Initialize count of pyrrole units
    pyrrole_units = 0

    # Iterate over rings
    for ring_atoms in ring_info.AtomRings():
        # Check if ring is size 5
        if len(ring_atoms) != 5:
            continue

        # Count nitrogen atoms in the ring
        num_nitrogen = 0
        for idx in ring_atoms:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 7:
                num_nitrogen += 1

        # If ring has exactly one nitrogen, count it as a pyrrole unit
        if num_nitrogen == 1:
            pyrrole_units += 1

    # Check if molecule is a polypyrrole
    if pyrrole_units >= 2:
        return True, f"Contains {pyrrole_units} pyrrole units"
    else:
        return False, f"Contains only {pyrrole_units} pyrrole unit(s). A polypyrrole requires two or more pyrrole units"


__metadata__ = {
    'chemical_class': {
        'name': 'polypyrrole',
        'definition': 'A compound composed of two or more pyrrole units.'
    }
}