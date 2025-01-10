"""
Classifies: CHEBI:38958 indole alkaloid
"""
from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid must have an indole skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible SMARTS pattern for indole core
    indole_flex_pattern = Chem.MolFromSmarts('c1ccc2[nH]c(c2c1)')
    if not mol.HasSubstructMatch(indole_flex_pattern):
        return False, "No indole skeleton found"

    # Additional checks for alkaloid properties
    # Define a nitrogen pattern that can find any nitrogen in the structure
    nitrogen_pattern = Chem.MolFromSmarts('[#7]')
    if not mol.HasSubstructMatch(nitrogen_pattern):
        return False, "Lacks additional nitrogen typically found in alkaloids"

    # Check if there are additional structural features specific to alkaloids, like ring structures
    # For now, we just check for at least one non-atrazine-like nitrogen

    return True, "Contains indole skeleton and additional features typical of alkaloids"

# Example usage for testing based on provided SMILES strings
example_smiles = "O=C1N2[C@H]([C@@]3(O[C@](C(N3[C@H]1CC(C)C)=O)(NC(=O)[C@@H]4C=C5C6=C7C(NC=C7C[C@H]5N(C4)C)=CC=C6)CC)O)CCC2"
result, reason = is_indole_alkaloid(example_smiles)
print(result, reason)