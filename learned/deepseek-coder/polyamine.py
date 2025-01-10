"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: CHEBI:88061 polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as any organic amino compound that contains two or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None  # Invalid SMILES string

    # Define amino group patterns: -NH2, -NH-, -NR2, [NH3+], [NH2+]
    amino_patterns = [
        Chem.MolFromSmarts("[NX3;H2]"),  # -NH2
        Chem.MolFromSmarts("[NX3;H1]"),  # -NH-
        Chem.MolFromSmarts("[NX3;H0]"),  # -NR2
        Chem.MolFromSmarts("[NH3+]"),    # [NH3+]
        Chem.MolFromSmarts("[NH2+]"),    # [NH2+]
    ]

    # Count unique amino groups
    amino_count = 0
    matched_atoms = set()  # To avoid counting the same atom multiple times
    for pattern in amino_patterns:
        if pattern is not None:  # Ensure the pattern is valid
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                atom_idx = match[0]  # Get the index of the nitrogen atom
                if atom_idx not in matched_atoms:
                    matched_atoms.add(atom_idx)
                    amino_count += 1

    # Check if there are at least two amino groups
    if amino_count >= 2:
        return True, f"Contains {amino_count} amino groups"
    else:
        return False, f"Contains only {amino_count} amino groups (need at least 2)"