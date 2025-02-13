"""
Classifies: CHEBI:15734 primary alcohol
"""
"""
Classifies: Primary Alcohol
"""
from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.
    A primary alcohol has a hydroxy group (-OH) attached to a carbon
    that is connected to either three hydrogen atoms or one other carbon and two hydrogen atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to get RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for primary alcohol
    # Identify carbon atoms bonded to OH group with only 1 other carbon (2 hydrogens)
    primary_alcohol_pattern_1 = Chem.MolFromSmarts("[C][OH]([H])[H]")
    # Optionally: a standalone primary alcohol carbon having three hydrogens
    primary_alcohol_pattern_2 = Chem.MolFromSmarts("[C][OH]([H])[H][H]")
    
    # Check for matches with the defined SMARTS patterns
    matches_1 = mol.GetSubstructMatches(primary_alcohol_pattern_1)
    matches_2 = mol.GetSubstructMatches(primary_alcohol_pattern_2)
    
    total_matches = len(matches_1) + len(matches_2)
    
    if total_matches > 0:
        return True, f"Contains {total_matches} primary alcohol group(s)."
    
    return False, "No primary alcohol group found."