"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone has an -OH group adjacent to a C=O group
    where the alpha carbon is a secondary carbon (attached to two carbons and has one hydrogen).

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

    # Define a new SMARTS pattern for secondary alpha-hydroxy ketones
    # The SMARTS pattern: a secondary carbon [C;R0X4] (having two carbon attachments) with an adjacent -OH and C=O
    pattern = Chem.MolFromSmarts("[#6]([C;R0X4;!R][OH])[C](=O)")
    
    # Check for match of the pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Contains a secondary alpha-hydroxy ketone functional group"
        
    return False, "Does not contain a secondary alpha-hydroxy ketone functional group"