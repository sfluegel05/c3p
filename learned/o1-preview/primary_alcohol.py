"""
Classifies: CHEBI:15734 primary alcohol
"""
"""
Classifies: CHEBI:15734 primary alcohol
"""

from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.
    A primary alcohol is a compound in which a hydroxy group (-OH) is attached to a saturated carbon atom
    which has either three hydrogen atoms attached to it (methanol), or only one other carbon atom and two hydrogen atoms attached to it (i.e., R-CH₂OH).
    
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
    
    # Define SMARTS pattern for primary alcohol
    # [#6;X4;H2,H3] - sp3 carbon with 2 or 3 hydrogens attached (primary carbon)
    # [OX2H] - hydroxyl group
    primary_alcohol_pattern = Chem.MolFromSmarts("[#6;X4;H2,H3][OX2H]")
    
    # Search for the primary alcohol pattern
    if mol.HasSubstructMatch(primary_alcohol_pattern):
        return True, "Contains primary alcohol group"
    else:
        return False, "Does not contain primary alcohol group"