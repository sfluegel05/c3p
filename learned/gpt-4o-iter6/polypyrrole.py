"""
Classifies: CHEBI:38077 polypyrrole
"""
from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is defined as a compound composed of two or more pyrrole units or
    pyrrole-like structures in larger, fused systems.

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

    # Define a list of SMARTS patterns for pyrrole and commonly recognized pyrrole-related structures
    pyrrole_patterns = [
        Chem.MolFromSmarts("n1cccc1"),    # Simple pyrrole ring
        Chem.MolFromSmarts("n1ccccc1"),   # Wider pyrrole-like fused ring
        Chem.MolFromSmarts("n3c4c5c2cnc6nc(c1c(ccc6)n12)c(c1c3c(c(c5)n4)c(cc1)n5)"),  # Porphyrin-like pattern
        Chem.MolFromSmarts("n2c1c(ncnc1ccc2)"),  # Azole-like pyrrole structure
    ]

    # Count the number of occurrences of pyrrole units in the molecule
    num_pyrrole_structures = 0

    for pattern in pyrrole_patterns:
        if pattern is None:
            continue

        matches = mol.GetSubstructMatches(pattern)
        num_pyrrole_structures += len(matches)

    # Check if there are two or more pyrrole/pyrrole-like units
    if num_pyrrole_structures >= 2:
        return True, f"Contains {num_pyrrole_structures} pyrrole-like structures"
    else:
        return False, "Contains less than two pyrrole-like structures"

# Example usage:
smiles_example = "C1=CNC=C1C1=CC=CN=C1"
is_polypyrrole(smiles_example)