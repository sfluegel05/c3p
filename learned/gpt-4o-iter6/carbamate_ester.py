"""
Classifies: CHEBI:23003 carbamate ester
"""
from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    Carbamate esters typically contain the carbamate functional group (-O-C(=O)-N-) or its derivatives.

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

    # Define SMARTS pattern for carbamate ester
    # This pattern matches standard carbamate esters
    carbamate_patterns = [
        Chem.MolFromSmarts("[#8]-[#6](=[#8])-[#7]")  # -O-C(=O)-N-
    ]

    # Check if molecule matches any carbamate ester pattern
    for pattern in carbamate_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains carbamate ester functional group"

    # No matching pattern found
    return False, "No carbamate ester functional group found"

# Test the function with different SMILES strings as needed.