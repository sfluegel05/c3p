"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is formed by the esterification of tetradecanoic acid with an alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define updated tetradecanoate ester functional group pattern
    # 13-carbon chain with ester linkage and additional carbon confirming alcohol linkage
    tetradecanoate_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCC(=O)O[CX4]")
    
    if mol.HasSubstructMatch(tetradecanoate_pattern):
        return True, "Contains tetradecanoate ester group"
    else:
        return False, "No tetradecanoate ester group found"