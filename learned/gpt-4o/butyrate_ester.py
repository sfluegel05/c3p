"""
Classifies: CHEBI:50477 butyrate ester
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester is defined as any carboxylic ester where the carboxylic acid component is butyric acid.

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

    # Define butyrate ester pattern (butyric acid part as ester: CCCC(=O)O)
    butyrate_pattern = Chem.MolFromSmarts("CCCC(=O)O")
    if not mol.HasSubstructMatch(butyrate_pattern):
        return False, "No butyrate ester group found"

    return True, "Contains a butyrate ester group"