"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: Alpha-Hydroxy Ketone
"""

from rdkit import Chem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is defined as a ketone containing a hydroxy group on 
    the alpha-carbon relative to the C=O group.

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
    
    # Define the alpha-hydroxy ketone SMARTS pattern:
    # [CX3](=O): Ketone carbon (sp2 hybridized carbon double-bonded to oxygen)
    # [CX4][OX2H]: Alpha-carbon (sp3 hybridized carbon) connected to a hydroxyl group
    alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[CX4][OX2H]")
    if alpha_hydroxy_ketone_pattern is None:
        return False, "Invalid SMARTS pattern"
    
    # Search for matches
    matches = mol.GetSubstructMatches(alpha_hydroxy_ketone_pattern)
    if matches:
        return True, "Contains an alpha-hydroxy ketone group"
    else:
        return False, "Does not contain an alpha-hydroxy ketone group"