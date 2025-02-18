"""
Classifies: CHEBI:15734 primary alcohol
"""
from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.

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

    # Look for a primary alcohol pattern: Hydroxyl group connected to a primary carbon
    # SMARTS pattern for primary alcohol: [CH2][OH]
    # This expression can be more complex to cover branching and other variations: [CH3][OH]
    primary_alcohol_pattern = Chem.MolFromSmarts("[CH2][OH]")

    # Check for matches of the primary alcohol pattern in the molecule
    if mol.HasSubstructMatch(primary_alcohol_pattern):
        return True, "Contains a primary alcohol group"

    # Also check for the case where the carbon is attached to three hydrogens and one hydroxyl: CH3OH
    methyl_alcohol_pattern = Chem.MolFromSmarts("[CH3][OH]")
    if mol.HasSubstructMatch(methyl_alcohol_pattern):
        return True, "Contains a methanol-like primary alcohol group"

    # If no primary alcohol pattern is found
    return False, "Does not contain a primary alcohol group"