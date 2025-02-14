"""
Classifies: CHEBI:23003 carbamate ester
"""
from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester is characterized by the functional group R¹–O–C(=O)–N(R²)–R³.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbamate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the carbamate ester SMARTS pattern
    carbamate_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)[OX2][CX4]")

    if carbamate_pattern is None:
        return False, "Invalid SMARTS pattern for carbamate ester"

    # Check if the molecule contains the carbamate ester pattern
    if mol.HasSubstructMatch(carbamate_pattern):
        return True, "Contains carbamate ester functional group"
    else:
        return False, "Does not contain carbamate ester functional group"