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
    
    # Define improved SMARTS pattern for alpha-hydroxy ketone
    # Ensure C-C(=O)-C-OH pattern is correctly identified
    alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts("[CX4][CX3](=O)[CX4][OX2H]")
    
    # Apply the pattern to the molecule
    if mol.HasSubstructMatch(alpha_hydroxy_ketone_pattern):
        return True, "Contains alpha-hydroxy ketone motif (C=O and adjacent C-OH)"
    
    # If pattern does not match, return false
    return False, "Does not contain alpha-hydroxy ketone motif"