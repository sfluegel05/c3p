"""
Classifies: CHEBI:46640 diketone
"""
from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone contains exactly two ketone functionalities.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diketone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a ketone functionality.
    # Ensure end atoms to prevent overlap bias in detection
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[CX3]")

    # Find all matches of ketone functionalities, ensuring proper atom constraints
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)

    # Improving the accuracy of ketone count by ensuring non-overlapping detection
    # Filtering the matches if there are any overlaps that might lead to double-counting
    unique_ketone_centers = set()
    for match in ketone_matches:
        unique_ketone_centers.add(match[1])  # Considering center of ketone (carbon)

    # Count the unique ketone functionalities by considering the core atom involved
    num_ketones = len(unique_ketone_centers)

    if num_ketones == 2:
        return True, "Molecule contains exactly two ketone functionalities"
    else:
        return False, f"Found {num_ketones} ketone functionalities, expected exactly 2"