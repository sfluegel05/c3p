"""
Classifies: CHEBI:17478 aldehyde
"""
from rdkit import Chem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.
    An aldehyde is defined as a compound containing a carbonyl group bonded to one hydrogen atom and one R group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the aldehyde group: R-C(=O)H
    aldehyde_pattern = Chem.MolFromSmarts("[#6]=O")

    # Check for the aldehyde group in the molecule
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Contains an aldehyde group (R-C(=O)H)"
    else:
        return False, "No aldehyde group found"

# Example usage
example_smiles = "C=O"  # Formyl group in SMILES
result, reason = is_aldehyde(example_smiles)
print(f"Is aldehyde: {result}, Reason: {reason}")