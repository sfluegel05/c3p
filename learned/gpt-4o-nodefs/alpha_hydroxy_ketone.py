"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
from rdkit import Chem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone has a hydroxyl group directly adjacent to a ketone group (C=O),
    i.e., R1-CH(OH)-C(=O)-R2 where R1 and R2 can be any group or hydrogen.

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
    
    # Define the alpha-hydroxy ketone SMARTS pattern (acyclic or cyclic)
    alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts("[C;H1](O)[C](=O)[C,c]")  # C-C(=O)-C with OH on first C

    # Check if the molecule has the alpha-hydroxy ketone substructure pattern
    if mol.HasSubstructMatch(alpha_hydroxy_ketone_pattern):
        return True, "Contains an alpha-hydroxy ketone group"
    else:
        return False, "Does not contain an alpha-hydroxy ketone group"