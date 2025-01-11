"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester contains the tetradecanoate moiety, which is a 14-carbon saturated fatty acid ester.

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

    # Look for tetradecanoate moiety pattern ("CCCCCCCCCCCCCC(=O)O")
    tetradecanoate_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)O")
    
    if mol.HasSubstructMatch(tetradecanoate_pattern):
        return True, "Contains tetradecanoate moiety"
    else:
        return False, "Does not contain tetradecanoate moiety"