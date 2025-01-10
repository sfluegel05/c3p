"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
from rdkit import Chem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is a ketone containing a hydroxy group on the alpha-carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more comprehensive SMARTS pattern for an alpha-hydroxy ketone
    # Pattern: carbonyl carbon, bonded to an alpha carbon with OH group
    alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts("[$([CX3](=O)[C-!!$(C=O,CO,N)])][$([CX4][OX2H])]")  # C(alpha)-C(=O)-C(alpha)-OH
    
    # Search for Kohen groups with adjacent OH
    substructures_found = mol.GetSubstructMatches(alpha_hydroxy_ketone_pattern)

    # Check if the pattern matches any substructure
    if substructures_found:
        return True, "Contains alpha-hydroxy ketone motif (C=O and adjacent C-OH)"
    
    # If pattern does not match, return false
    return False, "Does not contain alpha-hydroxy ketone motif"