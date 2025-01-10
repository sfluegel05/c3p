"""
Classifies: CHEBI:38958 indole alkaloid
"""
from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid must have an indole skeleton and typically contains other nitrogen atoms.

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

    # Define a flexible SMARTS pattern for the indole core
    indole_flex_pattern = Chem.MolFromSmarts('c1ccc2[nH]ccc2c1')  # Matches the indole core structure
    if not mol.HasSubstructMatch(indole_flex_pattern):
        return False, "No indole skeleton found"

    # Define pattern to capture any additional nitrogen atoms, which are typical in alkaloids
    additional_nitrogen_pattern = Chem.MolFromSmarts('[#7]')  # General pattern for nitrogen
    nitrogen_count = len(mol.GetSubstructMatches(additional_nitrogen_pattern))
    if nitrogen_count < 2:
        return False, f"Lacks sufficient nitrogen atoms for typical alkaloid structure, found {nitrogen_count}"

    return True, "Contains indole skeleton and features typical of alkaloids"

# Example usage for testing based on provided SMILES strings
example_smiles = "O=C1N2[C@H]([C@@]3(O[C@](C(N3[C@H]1CC(C)C)=O)(NC(=O)[C@@H]4C=C5C6=C7C(NC=C7C[C@H]5N(C4)C)=CC=C6)CC)O)CCC2"
result, reason = is_indole_alkaloid(example_smiles)
print(result, reason)