"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid contains exactly two carboxylic acid (COOH) groups.

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

    # Define the carboxylic acid pattern including extra constraints to avoid matches due to errors
    # Note: '?' allows for zero or one hydrogen in SMARTS to match OH/O/varied typo forms accurately.
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H0-,OX1H1]")

    # Get matches for the carboxylic acid pattern
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)

    # Count distinct groups, filtering to prevent gluings or internal errors/misidentifications in larger constructs
    distinct_matches = len(set(carboxylic_acid_matches))
    
    # Check the number of carboxylic acid groups
    if distinct_matches == 2:
        return True, "Contains exactly two carboxylic acid groups"
    elif distinct_matches > 2:
        return False, f"Contains {distinct_matches} carboxylic acid groups, which is more than 2"
    else:
        return False, f"Contains only {distinct_matches} carboxylic acid groups, less than or greater than needed"