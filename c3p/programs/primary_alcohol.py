"""
Classifies: CHEBI:15734 primary alcohol
"""
from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.
    A primary alcohol is defined as a hydroxyl group (-OH) attached to a terminal carbon,
    which is bonded to at least two hydrogen atoms or one carbon and two hydrogen atoms.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Enhance the pattern to capture primary alcohol connectivity including branching
    primary_alcohol_pattern = Chem.MolFromSmarts("[CX4;H2][OX2H]")
    
    # Check if the molecule has a substructure match for the enhanced primary alcohol pattern
    if mol.HasSubstructMatch(primary_alcohol_pattern):
        return True, "Contains a primary alcohol group"
    
    # Additional check for simple methyl alcohols as special cases
    methyl_alcohol_pattern = Chem.MolFromSmarts("[CX4;H3][OX2H]")
    if mol.HasSubstructMatch(methyl_alcohol_pattern):
        return True, "Contains a methanol-like primary alcohol group"
    
    # If no pattern matches, return false
    return False, "Does not contain a primary alcohol group"