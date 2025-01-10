"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid is any carboxylic acid containing exactly two carboxylic acid groups (either free or as part of different functional forms).

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
    
    # Define a broader pattern for carboxyl groups in various forms
    carboxyl_pattern_free = Chem.MolFromSmarts("C(=O)O")  # detects COOH (free) and COO-
    carboxyl_pattern_bound = Chem.MolFromSmarts("O=C([O,N])")  # broader match, catching amides or esters if necessary

    # Get matches for the broad patterns
    free_carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern_free)
    bound_carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern_bound)

    # Count the distinct set of carboxyl groups
    unique_carboxyl_set = set(free_carboxyl_matches + bound_carboxyl_matches)
    num_carboxyls = len(unique_carboxyl_set)

    # Dicarboxylic acids should have exactly two such groups
    if num_carboxyls == 2:
        return True, f"Contains exactly 2 carboxyl groups, indicating it is a dicarboxylic acid"
    elif num_carboxyls > 2:
        return False, f"Found {num_carboxyls} carboxyl groups; more than 2, cannot be a dicarboxylic acid"
    else:
        return False, f"Found {num_carboxyls} carboxyl groups; need exactly 2"

# The broad pattern approach ensures we account for free COOH, COO- groups, and those bound in other chemical functionalities.