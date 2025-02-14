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
    A primary alcohol has a hydroxy group (-OH) attached to a carbon with
    either three hydrogen atoms or bonded to one carbon and two hydrogen atoms.

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
    
    # Define SMARTS pattern for primary alcohols
    # Primary carbon with OH pertaining to primary alcohols
    primary_alcohol_pattern = Chem.MolFromSmarts("[CX4;!R;H2,H3][OH]")

    # Check for matches with the defined SMARTS pattern
    matches = mol.GetSubstructMatches(primary_alcohol_pattern)
    
    if matches:
        return True, f"Contains {len(matches)} primary alcohol group(s)."

    return False, "No primary alcohol group found."