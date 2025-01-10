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
    
    # Define the SMARTS pattern for an alpha-hydroxy ketone:
    # - Recognize a carbonyl group "[CX3](=O)"
    # - Check for an adjacent carbon (alpha-carbon) bonded to an OH group "[CX4][CX4][OX2H]"
    alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[CX4H1][OX2H]")

    # Ensure the SMARTS pattern is valid
    if alpha_hydroxy_ketone_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Search for the alpha-hydroxy ketone pattern in the molecule
    substructures_found = mol.GetSubstructMatches(alpha_hydroxy_ketone_pattern)

    # Check if the pattern matches any substructure
    if substructures_found:
        return True, "Contains alpha-hydroxy ketone motif (C=O and adjacent C-OH)"
    
    # If pattern does not match, return false
    return False, "Does not contain alpha-hydroxy ketone motif"