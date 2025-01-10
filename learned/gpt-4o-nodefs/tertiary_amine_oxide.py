"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
from rdkit import Chem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    A tertiary amine oxide is characterized by a nitrogen atom bonded to three carbon atoms
    and an oxygen atom (N+ → O-).

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

    # Define tertiary amine oxide pattern: N(+)-O(-)
    amine_oxide_pattern = Chem.MolFromSmarts("[N+](C)(C)(C)[O-]")

    # Check for the pattern in the molecule
    if not mol.HasSubstructMatch(amine_oxide_pattern):
        return False, "No tertiary amine oxide structure found"

    return True, "Contains a tertiary amine oxide structure (N+ → O-)"

# Test the function with an example
smiles_example = "C[N+](C)([O-])C"
is_tertiary_amine_oxide(smiles_example)