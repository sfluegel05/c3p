"""
Classifies: CHEBI:46640 diketone
"""
from rdkit import Chem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone is a compound that contains two ketone functionalities.

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

    # Define ketone SMARTS pattern: carbonyl carbon double-bonded to oxygen and single-bonded to two carbons
    ketone_pattern = Chem.MolFromSmarts('[#6;X3](=O)([#6])[#6]')
    if ketone_pattern is None:
        return False, "Invalid ketone pattern"

    # Find ketone matches
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    num_ketones = len(ketone_matches)

    if num_ketones == 2:
        return True, "Contains exactly two ketone functionalities"
    else:
        return False, f"Contains {num_ketones} ketone functionalities, need exactly 2"