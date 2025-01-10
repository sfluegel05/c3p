"""
Classifies: CHEBI:24128 furanocoumarin
"""
from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin consists of a furan ring fused with a coumarin, allowing for different fusion patterns.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for different types of furanocoumarin structures
    linear_furanocoumarin_pattern = Chem.MolFromSmarts('c1co2c(=O)ccc2cc1')  # Generic pattern for linear furanocoumarins
    angular_furanocoumarin_pattern = Chem.MolFromSmarts('c1cc2coc(=O)c2cc1')  # Generic pattern for angular furanocoumarins
    
    # Check for the presence of any recognized furanocoumarin core structure
    if mol.HasSubstructMatch(linear_furanocoumarin_pattern):
        return True, "Contains furanocoumarin core (linear pattern)"
    if mol.HasSubstructMatch(angular_furanocoumarin_pattern):
        return True, "Contains furanocoumarin core (angular pattern)"

    return False, "No furanocoumarin core structure found"