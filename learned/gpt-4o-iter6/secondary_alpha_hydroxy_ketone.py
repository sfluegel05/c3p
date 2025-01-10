"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone has an -OH group adjacent to a C=O group
    where the central carbon is secondary (has one hydrogen and one other organic group).

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

    # Look for a secondary alpha-hydroxy ketone pattern,
    # specifically, a secondary carbon linked with an OH group and adjacent to a C=O group
    # SMARTS pattern corrections for capturing a secondary alpha-hydroxy ketone structure
    pattern = Chem.MolFromSmarts("[#6;R1;$([CH1]([C,c])[OX2H])[C](=O)]")
    
    # Check for match of the pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Contains a secondary alpha-hydroxy ketone functional group"
        
    return False, "Does not contain a secondary alpha-hydroxy ketone functional group"