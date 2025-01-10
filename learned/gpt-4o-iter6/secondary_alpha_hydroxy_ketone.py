"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone has an -OH group adjacent to a C=O group
    where the alpha carbon is a secondary carbon (has one hydrogen and at least
    one organyl group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for secondary alpha-hydroxy ketones
    # A secondary carbon with an OH group and adjacent to a carbonyl group (C(=O))
    # Here, we look for [C] indicating a carbon atom with one hydrogen and at least one attachment
    pattern = Chem.MolFromSmarts("[#6][CH]([#6])[OH][C](=O)") 
    
    # Check for match of the pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Contains a secondary alpha-hydroxy ketone functional group"
        
    return False, "Does not contain a secondary alpha-hydroxy ketone functional group"