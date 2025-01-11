"""
Classifies: CHEBI:16385 organic sulfide
"""
from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide based on its SMILES string.
    An organic sulfide is characterized by a sulfur atom bonded to two carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sulfur atom bonded to two carbon atoms (basic thioether structure)
    sulfide_pattern = Chem.MolFromSmarts("[CX4][SX2][CX4]")  # Simple pattern for sulfide
    if mol.HasSubstructMatch(sulfide_pattern):
        return True, "Contains sulfur atom bonded to two carbon atoms"

    # Check for more complex organic sulfide structures
    complex_sulfide_pattern = Chem.MolFromSmarts("[SX2]")  # Matches sulfur with 2 single bonds
    if mol.HasSubstructMatch(complex_sulfide_pattern):
        return True, "Contains a sulfur atom likely to be part of a complex organic sulfide"

    return False, "Does not contain sulfur bound to two carbon atoms typical of organic sulfides"