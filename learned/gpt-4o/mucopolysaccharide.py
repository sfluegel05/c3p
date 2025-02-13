"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is likely to be a mucopolysaccharide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule has characteristics of a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define uronic acid patterns
    uronic_acid_pattern = Chem.MolFromSmarts("[C;R0](=O)[O;R0][CH][CH]")  # Carboxyl group on non-ring C adjacent to sugars
    # Define glycosamine patterns
    glycosamine_pattern = Chem.MolFromSmarts("[CH][CH][NX3][CH2]")  # Nitrogen attached to sugar-like structure
    
    # Define sulfate esterification pattern
    sulfate_pattern = Chem.MolFromSmarts("[#6][O,SX4](=O)(=O)[O,S][#6]")  # Consider sulfate attached to sugars

    # Check for the presence of patterns
    uronic_acid_matches = mol.HasSubstructMatch(uronic_acid_pattern)
    glycosamine_matches = mol.HasSubstructMatch(glycosamine_pattern)
    sulfate_matches = mol.HasSubstructMatch(sulfate_pattern)

    if uronic_acid_matches and glycosamine_matches:
        if sulfate_matches:
            return True, "Characteristics of mucopolysaccharide, including sulfate esterification"
        else:
            return True, "Characteristics of mucopolysaccharide without detectable sulfate groups"
    
    return False, "Missing key features of a mucopolysaccharide"