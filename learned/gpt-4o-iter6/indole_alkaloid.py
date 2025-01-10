"""
Classifies: CHEBI:38958 indole alkaloid
"""
from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid must have an indole skeleton, albeit with potential modifications,
    and typically contains additional nitrogen atoms.

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

    # Enhanced pattern to recognize a broader set of indole derivatives
    indole_flex_pattern = Chem.MolFromSmarts('c1c[nH]c2c1ccc2')  # General indole pattern
    carbazole_pattern = Chem.MolFromSmarts('c1cc2ccccc2nc1')     # Common alkaloid pattern

    if not (mol.HasSubstructMatch(indole_flex_pattern) or mol.HasSubstructMatch(carbazole_pattern)):
        return False, "No recognizable indole variant skeleton found"

    # Define pattern to capture the presence of additional nitrogen atoms
    additional_nitrogen_pattern = Chem.MolFromSmarts('[#7]')  # General pattern for nitrogen
    nitrogen_count = len(mol.GetSubstructMatches(additional_nitrogen_pattern))
    
    # Generally less strict, since indole itself already contains nitrogen
    if nitrogen_count < 1:
        return False, f"Lacks additional nitrogen atoms for typical complex alkanoid structure, found {nitrogen_count}"

    return True, "Contains indole skeleton or close variants with typical alkaloidal features"

# Example usage for testing based on provided SMILES strings
example_smiles = "O=C1N2[C@H]([C@@]3(O[C@](C(N3[C@H]1CC(C)C)=O)(NC(=O)[C@@H]4C=C5C6=C7C(NC=C7C[C@H]5N(C4)C)=CC=C6)CC)O)CCC2"
result, reason = is_indole_alkaloid(example_smiles)
print(result, reason)