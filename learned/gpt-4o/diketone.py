"""
Classifies: CHEBI:46640 diketone
"""
from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone contains exactly two ketone (C=O) functionalities.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for a ketone group
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[CX4]")

    # Find ketone groups in the molecule
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    # Check if there are exactly two ketone groups
    if len(ketone_matches) == 2:
        return True, "Contains exactly two ketone groups, hence a diketone."
    else:
        return False, f"Contains {len(ketone_matches)} ketone groups, not a diketone."

# Examples
# Test the function with a diketone example SMILES
result, reason = is_diketone("CCCCCC(=O)CC(=O)CCCCCC")
print(result, reason)

# Test the function with a non-diketone example SMILES
result, reason = is_diketone("CCCCCCC=O")
print(result, reason)