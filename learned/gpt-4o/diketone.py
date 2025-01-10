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

    # Define the SMARTS pattern for a ketone functionality
    ketone_pattern = Chem.MolFromSmarts("C(=O)C")

    # Find all matches of ketone functionalities
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)

    # Count the number of ketone functionalities found
    num_ketones = len(ketone_matches)

    if num_ketones == 2:
        return True, "Molecule contains exactly two ketone functionalities"
    else:
        return False, f"Found {num_ketones} ketone functionalities, expected exactly 2"