"""
Classifies: CHEBI:46640 diketone
"""
from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone is defined as a compound that contains two ketone functionalities.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is considered a diketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ketone group pattern: a carbonyl group (C=O) with a carbon atom on at least one side
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[#6]")
    if ketone_pattern is None:
        return False, "Invalid ketone SMARTS pattern"

    # Find substructure matches for ketone groups
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)

    # Check for exactly two ketone groups, considering all detected ketones
    if len(ketone_matches) == 2:
        return True, f"Contains exactly 2 ketone groups, sufficient for diketone classification"

    return False, f"Found {len(ketone_matches)} ketone groups, exactly 2 needed for diketone"

# Example usage:
# print(is_diketone("CC(=O)CC(C)=O"))  # acetylacetone: SMILES: CC(=O)CC(C)=O