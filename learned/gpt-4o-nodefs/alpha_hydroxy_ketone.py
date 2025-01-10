"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
from rdkit import Chem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone has a hydroxyl group adjacent to a ketone group (C=O).

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
    
    # Define alpha-hydroxy ketone pattern
    alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts("O[C;R0][C;R0](=O)")
    
    # Check if the molecule has the alpha-hydroxy ketone substructure
    if mol.HasSubstructMatch(alpha_hydroxy_ketone_pattern):
        return True, "Contains an alpha-hydroxy ketone group"
    else:
        return False, "Does not contain an alpha-hydroxy ketone group"

# Example usage:
# smiles = "CC(=O)CO" # hydroxyacetone
# is_alpha_hydroxy_ketone(smiles)