"""
Classifies: CHEBI:23003 carbamate ester
"""
from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester is an ester of carbamic acid (\(R^1OC(=O)NR^2R^3\)).

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

    # SMARTS pattern for carbamate ester: O-C(=O)-N
    carbamate_pattern = Chem.MolFromSmarts("O=C(O)N")
    
    # Check for carbamate ester substructure
    if mol.HasSubstructMatch(carbamate_pattern):
        return True, "Contains the carbamate ester group"

    return False, "Does not contain a carbamate ester group"