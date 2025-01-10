"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is formed by the esterification of tetradecanoic acid (C13H27COOH).

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
    
    # Define SMARTS pattern for tetradecanoic acid esterified group: C13H27C(=O)O-
    tetradecanoate_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCC(=O)O")
    
    # Check for presence of tetradecanoate ester functionality
    if not mol.HasSubstructMatch(tetradecanoate_pattern):
        return False, "No tetradecanoate ester group found"
    
    # Verify if multiple such groups are present
    ester_matches = mol.GetSubstructMatches(tetradecanoate_pattern)
    if len(ester_matches) < 1:
        return False, "Missing tetradecanoate esters, none detected"

    return True, "Contains tetradecanoate ester functionality"