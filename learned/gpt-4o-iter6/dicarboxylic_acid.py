"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid is any carboxylic acid containing exactly two carboxyl groups.

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
    
    # Define pattern for free carboxyl group as -C(=O)O
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if carboxyl_pattern is None:
        return False, "Invalid carboxyl pattern"

    # Get matches for carboxyl groups
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)

    # Count the distinct carboxyl groups
    num_carboxyl_groups = len(carboxyl_matches)

    # Check for exactly two carboxyl groups for classification
    if num_carboxyl_groups == 2:
        return True, "Contains exactly 2 carboxyl groups, indicating it's a dicarboxylic acid"
    elif num_carboxyl_groups < 2:
        return False, f"Found {num_carboxyl_groups} carboxyl groups; need exactly 2"
    else:
        return False, f"Found {num_carboxyl_groups} carboxyl groups; more than 2, cannot be a dicarboxylic acid"

# By focusing on distinct carboxyl groups that fit the accurate description of isolated acid functionalities, the classification should better line up with known data.