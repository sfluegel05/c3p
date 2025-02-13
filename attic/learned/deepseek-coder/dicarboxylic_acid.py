"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: CHEBI:17855 dicarboxylic acid
"""
from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid contains exactly two carboxylic acid groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dicarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    
    # Find all matches of the carboxylic acid pattern
    matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    
    # Check if there are exactly 2 carboxylic acid groups
    if len(matches) == 2:
        return True, "Contains exactly two carboxylic acid groups"
    else:
        return False, f"Found {len(matches)} carboxylic acid groups, need exactly 2"