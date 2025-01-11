"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol is a compound that contains exactly two alcoholic hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for alcoholic hydroxyl groups (-OH attached to a carbon)
    alcoholic_hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    
    # Find all matches for the alcoholic hydroxyl pattern
    alcoholic_hydroxyl_matches = mol.GetSubstructMatches(alcoholic_hydroxyl_pattern)
    
    # Count the number of alcoholic hydroxyl groups
    alcoholic_hydroxyl_count = len(alcoholic_hydroxyl_matches)

    # Check if the molecule has exactly two alcoholic hydroxyl groups
    if alcoholic_hydroxyl_count == 2:
        return True, "Contains exactly two alcoholic hydroxyl groups"
    else:
        return False, f"Contains {alcoholic_hydroxyl_count} alcoholic hydroxyl groups, need exactly 2"