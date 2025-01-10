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
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for alpha-hydroxy ketone pattern
    # Carbon (not necessarily in a ring) with C=O adjacent to C-OH
    alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts("[#6][C](=O)[C][O]")
    if mol.HasSubstructMatch(alpha_hydroxy_ketone_pattern):
        return True, "Contains alpha-hydroxy ketone motif (C=O and adjacent C-OH)"
    
    return False, "Does not contain alpha-hydroxy ketone motif"