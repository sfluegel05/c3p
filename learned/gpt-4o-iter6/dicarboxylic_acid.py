"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid is any carboxylic acid containing exactly two free carboxyl groups.

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

    # Define a pattern for free carboxyl groups (COOH)
    carboxyl_pattern_free = Chem.MolFromSmarts("C(=O)[OH]")
    free_carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern_free)
    
    # Count number of free carboxyl groups
    num_free_carboxyls = len(free_carboxyl_matches)

    # Dicarboxylic acids should have exactly two free carboxyl groups
    if num_free_carboxyls == 2:
        return True, f"Contains exactly 2 free carboxyl groups, indicating it is a dicarboxylic acid"
    elif num_free_carboxyls > 2:
        return False, f"Found {num_free_carboxyls} free carboxyl groups; more than 2, cannot be a dicarboxylic acid"
    else:
        return False, f"Found {num_free_carboxyls} free carboxyl groups; need exactly 2"

# Note: This implementation assumes the interest in identifying free carboxyl groups
# avoiding esters or amides. Additional accuracy can require structural context or 
# experimental validation in complex systems.