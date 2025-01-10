"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
from rdkit import Chem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    A tertiary amine oxide is characterized by a nitrogen atom bonded to three carbon atoms and an oxygen atom (N+ → O-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine oxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define more comprehensive patterns for tertiary amine oxides
    patterns = [
        Chem.MolFromSmarts("[N+](C)(C)(C)[O-]"), # Basic tertiary amine oxide
        Chem.MolFromSmarts("[N+](C)(C)(C)O"), # Allow O without formal charge
        Chem.MolFromSmarts("[N+]([R][C])(C)[O-]"), # Ring-substituted
        Chem.MolFromSmarts("[N+](C)[C](C)[O-]"), # Close proximity of O
        Chem.MolFromSmarts("[N+](C)([O-])[R]"), # Alternative ring positions
        Chem.MolFromSmarts("[R][N+]([R])([R])[O-]"), # Extensive ring systems
        Chem.MolFromSmarts("[N+](C)(C)[O-]C"), # Variant of configurations
        Chem.MolFromSmarts("[C][N+]([C])([C])[O-]R") # More general contexts
    ]

    # Check for any of the patterns in the molecule
    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a variant of the tertiary amine oxide structure (N+ → O-)"

    return False, "No tertiary amine oxide structure found"

# Test the function with an example
smiles_example = "C[N+](C)([O-])C"
print(is_tertiary_amine_oxide(smiles_example))