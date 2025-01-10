"""
Classifies: CHEBI:50477 butyrate ester
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester contains a butyric acid moiety attached via an ester linkage.

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

    # Look for butyrate ester pattern (CCCC(=O)O)
    butyrate_ester_pattern = Chem.MolFromSmarts("CCCC(=O)O")
    if not mol.HasSubstructMatch(butyrate_ester_pattern):
        return False, "No butyrate ester functional group found"

    return True, "Contains butyrate ester functional group"