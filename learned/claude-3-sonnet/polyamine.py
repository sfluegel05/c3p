"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: CHEBI:26144 polyamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is an organic compound containing two or more amino groups.

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

    # Check if molecule is organic (contains carbon)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Not an organic compound"

    # Define patterns for different types of nitrogen groups
    patterns = {
        # Valid amine patterns
        'primary_amine': '[NX3H2][CX4]',  # Aliphatic primary amine
        'secondary_amine': '[NX3H1]([CX4])[CX4]',  # Aliphatic secondary amine
        'tertiary_amine': '[NX3H0]([CX4])([CX4])[CX4]',  # Aliphatic tertiary amine
        'aromatic_amine': '[NX3H2]c',  # Aromatic primary amine
        'protonated_amine': '[NX4+]',  # Protonated amine
        
        # Patterns to exclude
        'amide': '[NX3][CX3](=[OX1])',  # Amide
        'nitro': '[$([NX3](=O)=O),$([NX3+](=O)[O-])]',  # Nitro
        'aromatic_N': '[nX2,nX3]',  # Aromatic nitrogen in ring
        'guanidine': '[NX3][CX3](=[NX2])[NX3]',  # Guanidine group
        'imine': '[NX2]=[CX3]',  # Imine
    }

    # Convert patterns to RDKit molecules
    patterns = {name: Chem.MolFromSmarts(pattern) for name, pattern in patterns.items()}

    # Count different types of amino groups
    counts = {}
    for name, pattern in patterns.items():
        counts[name] = len(mol.GetSubstructMatches(pattern))

    # Calculate total valid amine groups
    total_amines = (
        counts['primary_amine'] +
        counts['secondary_amine'] +
        counts['tertiary_amine'] +
        counts['aromatic_amine'] +
        counts['protonated_amine']
    )

    # Special handling for molecules with both amides and amines
    if total_amines >= 2:
        return True, (f"Contains {total_amines} amino groups "
                     f"(Primary: {counts['primary_amine'] + counts['aromatic_amine']}, "
                     f"Secondary: {counts['secondary_amine']}, "
                     f"Tertiary: {counts['tertiary_amine']}, "
                     f"Protonated: {counts['protonated_amine']})")

    # If we have an amide and at least one other amine group, count it
    elif total_amines == 1 and counts['amide'] >= 1:
        return True, (f"Contains 2 amino groups (including 1 amide nitrogen)")

    return False, f"Contains only {total_amines} amino groups, minimum 2 required"