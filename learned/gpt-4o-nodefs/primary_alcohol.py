"""
Classifies: CHEBI:15734 primary alcohol
"""
from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.
    A primary alcohol has the hydroxyl group (-OH) connected to a primary carbon.

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
    
    # Primary alcohol pattern: carbon with hydroxyl group attached and only one other carbon connected
    # [CH2][OH] is the pattern indicating a primary alcohol; CH2 or CH3-CH(OH)R where R is not a carbon atom
    primary_alcohol_pattern = Chem.MolFromSmarts("[CH2][OH]")
    
    # Check for matches of the primary alcohol pattern
    if mol.HasSubstructMatch(primary_alcohol_pattern):
        return True, "Structure contains a primary alcohol group"
    
    return False, "No primary alcohol group found"