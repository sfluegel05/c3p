"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: CHEBI:2468 alpha-hydroxy ketone
"""
from rdkit import Chem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is a ketone with a hydroxyl group on the alpha-carbon relative to the C=O group.

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

    # Define the SMARTS pattern for an alpha-hydroxy ketone
    # The pattern looks for a carbonyl group (C=O) with a hydroxyl group (OH) on the alpha-carbon
    alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[CX4][OX2H]")

    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(alpha_hydroxy_ketone_pattern):
        return True, "Contains a ketone with a hydroxyl group on the alpha-carbon"
    else:
        return False, "No alpha-hydroxy ketone pattern found"