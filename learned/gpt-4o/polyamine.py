"""
Classifies: CHEBI:88061 polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine contains two or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count amino groups (nitrogen atoms, N)
    amine_count = 0
    for atom in mol.GetAtoms():
        # Check if the atom is nitrogen with two or fewer other atoms connected
        if atom.GetAtomicNum() == 7:  # Nitrogen has atomic number 7
            # Count the number of connected non-hydrogen atoms
            # We want to ensure it is primary/secondary amine type
            non_h_count = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() != 1)
            if non_h_count <= 2:
                amine_count += 1

    # Determine classification
    if amine_count >= 2:
        return True, f"Molecule has {amine_count} amino groups, classifying as polyamine"
    else:
        return False, f"Molecule has {amine_count} amino group(s), not enough to classify as polyamine"