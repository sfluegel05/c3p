"""
Classifies: CHEBI:23003 carbamate ester
"""
from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester is defined by the presence of the carbamate functional group (-O-C(=O)-N-).

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

    # Look for carbamate ester group pattern (-O-C(=O)-N-)
    carbamate_pattern = Chem.MolFromSmarts("[$([NX3]([#6])[#6])][CX3](=[OX1])[OX2]")
    
    if mol.HasSubstructMatch(carbamate_pattern):
        return True, "Contains carbamate ester functional group"

    return False, "No carbamate ester functional group found"