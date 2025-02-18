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
    
    # Define SMARTS patterns for primary, secondary, and tertiary amines
    primary_amine_pattern = Chem.MolFromSmarts("[NX3;H2][C]")
    secondary_amine_pattern = Chem.MolFromSmarts("[NX3;H1][C][C]")
    tertiary_amine_pattern = Chem.MolFromSmarts("[NX3][C][C][C]")
    
    # Match SMARTS patterns
    if mol.HasSubstructMatch(primary_amine_pattern):
        return True, "Contains primary amine group"
    if mol.HasSubstructMatch(secondary_amine_pattern):
        return True, "Contains secondary amine group"
    if mol.HasSubstructMatch(tertiary_amine_pattern):
        return True, "Contains tertiary amine group"
    
    return False, "No amine group found, or nitrogen not bonded as required for amines"