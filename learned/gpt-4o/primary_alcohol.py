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
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for primary alcohol
    # [CX4H2][OH] means a carbon with two hydrogens attached to an OH group
    primary_alcohol_pattern = Chem.MolFromSmarts("[CX4;H2][OH]")

    # Check for primary alcohol match
    matches = mol.GetSubstructMatches(primary_alcohol_pattern)
    if matches:
        return True, f"Contains {len(matches)} primary alcohol group(s)."
    
    return False, "No primary alcohol group found."