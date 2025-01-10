"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
from rdkit import Chem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    A tertiary amine oxide is characterized by a nitrogen atom bonded to three carbon atoms
    and an oxygen atom (N+ → O-), possibly within a diverse array of structural contexts.

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

    # Define tertiary amine oxide patterns
    # Pattern for tertiary amine oxide with various structural contexts
    patterns = [
        Chem.MolFromSmarts("[N+](C)(C)(C)[O-]"), # Exact match
        Chem.MolFromSmarts("[N+](C)(C)(C)O"), # Allow unconjugated O
        Chem.MolFromSmarts("[N+](C)(C)C[O-]"), # Oxygen bonded to outer carbon
        Chem.MolFromSmarts("[R][N+](C)(C)[O-]"), # Incorporate ring
        Chem.MolFromSmarts("[N+](C)(C)[O-]C"), # Different configurations with O
    ]

    # Check for the pattern in the molecule
    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a variant of the tertiary amine oxide structure (N+ → O-)"

    return False, "No tertiary amine oxide structure found"

# Test the function with an example
smiles_example = "C[N+](C)([O-])C"
print(is_tertiary_amine_oxide(smiles_example))