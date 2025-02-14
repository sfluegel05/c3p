"""
Classifies: CHEBI:23003 carbamate ester
"""
from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester is defined as any ester of carbamic acid or its N-substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbamate ester, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define compact SMARTS pattern for identifying carbamate ester
    # Pattern: O=C(N)OC, generalized to account for substitutions
    carbamate_ester_pattern = Chem.MolFromSmarts("O=C(N[*])OC[*]")

    # Check for the carbamate ester presence using the defined pattern
    if mol.HasSubstructMatch(carbamate_ester_pattern):
        return True, "Contains carbamate ester functional group"

    return False, "Does not contain carbamate ester functional group"