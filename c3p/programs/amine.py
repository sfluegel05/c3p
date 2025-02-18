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
    # Nitrogen atom with at least one hydrocarbyl group attached.
    amine_pattern = Chem.MolFromSmarts("[NX3;!$(NC=[O,S,N])][C]")
    
    # Match SMARTS patterns
    if mol.HasSubstructMatch(amine_pattern):
        return True, "Contains amine group"
    
    return False, "No amine group found, or nitrogen not bonded as required for amines"