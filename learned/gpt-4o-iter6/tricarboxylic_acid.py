"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: tricarboxylic acid
"""
from rdkit import Chem

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid is defined as having three carboxylic acid groups (-COOH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define carboxylic acid group using SMARTS
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    
    # Count the number of carboxylic acid groups
    n_carboxylic_acids = len(carboxylic_acid_matches)
    if n_carboxylic_acids == 3:
        return True, "Contains exactly three carboxylic acid groups"
    else:
        return False, f"Contains {n_carboxylic_acids} carboxylic acid groups, expected exactly 3"