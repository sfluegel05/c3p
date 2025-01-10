"""
Classifies: CHEBI:23003 carbamate ester
"""
from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester is characterized by the presence of the carbamate functional group (-O-C(=O)-N-).

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

    # Define SMARTS pattern for carbamate ester: Includes nitrogen both in and out of rings
    # This checks for the -O-C(=O)-N- group and allows variation in surrounding structure
    carbamate_pattern = Chem.MolFromSmarts("[O]-C(=O)-[N]")

    if mol.HasSubstructMatch(carbamate_pattern):
        return True, "Contains carbamate ester functional group"

    # No matching pattern found
    return False, "No carbamate ester functional group found"