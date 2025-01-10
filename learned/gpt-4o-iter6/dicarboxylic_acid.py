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

    # Carboxylic acid pattern (distinct free carboxyl groups only)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)

    # Ensure the carboxyl groups are not part of esters or amides 
    # by checking connection(s) specifically at the 'O' of the carboxyl group
    carboxyl_free_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    free_matches = mol.GetSubstructMatches(carboxyl_free_pattern)
    
    # Count the distinct carboxylic acid groups that aren't bound in esters/amides
    num_true_carboxyls = len(free_matches)
    
    # Expect exactly 2 individual free carboxyl groups, not part of esters or amides
    if num_true_carboxyls == 2:
        return True, "Contains exactly 2 free carboxyl groups, indicating it is a dicarboxylic acid"
    elif num_true_carboxyls > 2:
        return False, f"Found {num_true_carboxyls} carboxyl groups; more than 2, cannot be a dicarboxylic acid"
    else:
        return False, f"Found {num_true_carboxyls} carboxyl groups; need exactly 2"

# Note: While the code aims to be comprehensive, real-case implementation should include 
# comprehensive molecular error handling, broader testing, and further feature specifics as 
# often these entities are part of larger molecular architectures, not single forms.