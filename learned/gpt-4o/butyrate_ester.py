"""
Classifies: CHEBI:50477 butyrate ester
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester contains a butyric acid component (CCCC(=O)O) esterified with an alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Butyric acid ester linkage pattern, allowing variable groups after the ester oxygen [C;H2,C;H3,C#C,...]
    # We want to ensure "CCCC(=O)O" is part of an ester
    butyrate_ester_pattern = Chem.MolFromSmarts("CCCC(=O)OC")
    
    # Perform the match
    if mol.HasSubstructMatch(butyrate_ester_pattern):
        return True, "Contains butyric acid ester linkage"

    return False, "No butyric acid ester linkage found"