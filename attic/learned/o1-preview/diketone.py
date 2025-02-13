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
        bool: True if the molecule is a diketone, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for ketone functional group
    # Carbonyl carbon (C=O) where the carbon is single-bonded to any carbon
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[#6]")
    if ketone_pattern is None:
        return False, "Invalid ketone SMARTS pattern"

    # Find all ketone matches in the molecule
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    num_ketones = len(ketone_matches)

    if num_ketones >= 2:
        return True, f"Contains {num_ketones} ketone functionalities"
    else:
        return False, f"Contains {num_ketones} ketone groups, need at least 2 for diketone"