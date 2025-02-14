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

    # Corrected carbamate ester SMARTS pattern (O-C(=O)-N)
    carbamate_smarts = "[OX2H0][CX3](=O)[NX3;!$(N-C=O)]"

    # Create the SMARTS pattern molecule
    carbamate_pattern = Chem.MolFromSmarts(carbamate_smarts)
    if carbamate_pattern is None:
        return False, "Invalid SMARTS pattern for carbamate ester"

    # Check if the molecule contains the carbamate ester pattern
    if mol.HasSubstructMatch(carbamate_pattern):
        return True, "Contains carbamate ester functional group"
    else:
        return False, "Does not contain carbamate ester functional group"