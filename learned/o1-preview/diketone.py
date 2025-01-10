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
        bool: True if molecule is a diketone, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the ketone functional group pattern
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[CX4H3,CX3H2,CX3H,CX3]")  # Carbonyl carbon double-bonded to oxygen and single-bonded to carbon
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    num_ketones = len(ketone_matches)

    if num_ketones == 2:
        return True, "Contains exactly two ketone functionalities"
    else:
        return False, f"Contains {num_ketones} ketone groups, need exactly 2 for diketone"