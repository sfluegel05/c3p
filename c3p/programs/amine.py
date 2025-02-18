"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine is a compound formally derived from ammonia by replacing one,
    two, or three hydrogen atoms by hydrocarbyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Refined SMARTS patterns for primary, secondary, and tertiary amines:
    # Exclude amines that are part of other functionalities like amide, nitro, and other complex bonds
    primary_amine_pattern = Chem.MolFromSmarts("[NX3H2][C]")
    secondary_amine_pattern = Chem.MolFromSmarts("[NX3H][C][C]")
    tertiary_amine_pattern = Chem.MolFromSmarts("[NX3]([C])[C][C]")
    
    # Combine the patterns and check the molecule
    if (mol.HasSubstructMatch(primary_amine_pattern) or
        mol.HasSubstructMatch(secondary_amine_pattern) or
        mol.HasSubstructMatch(tertiary_amine_pattern)):
        return True, "Contains amine group"
    
    # Check if nitrogen is part of excluded functionality
    excluded_functionalities = Chem.MolFromSmarts("[NX3;R]")  # Example for framework rings, can add more
    if mol.HasSubstructMatch(excluded_functionalities):
        return False, "Nitrogen part of non-amine functional group"
    
    return False, "No amine group found, or nitrogen not bonded as required for amines"