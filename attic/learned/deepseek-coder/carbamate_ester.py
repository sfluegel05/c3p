"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies: CHEBI:23004 carbamate ester
"""
from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester is any ester of carbamic acid or its N-substituted derivatives.

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

    # Define the carbamate ester pattern: -O-C(=O)-N-
    carbamate_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[NX3]")
    
    # Check if the pattern exists in the molecule
    if mol.HasSubstructMatch(carbamate_pattern):
        return True, "Contains the carbamate ester functional group (-O-C(=O)-N-)"
    else:
        return False, "Does not contain the carbamate ester functional group (-O-C(=O)-N-)"