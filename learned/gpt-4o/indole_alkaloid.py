"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: Indole Alkaloid
"""
from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid must contain an indole skeleton (a bicyclic structure consisting
    of a benzene ring fused to a pyrrole ring) and also contain functional groups
    or structural features characteristic of alkaloids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general indole skeleton pattern with a broader nitrogen definition
    indole_pattern = Chem.MolFromSmarts('c1ccc2c(c1)cnc2')  # A simple fused bicyclic structure

    # General patterns for additional nitrogen-containing functionalities typical in alkaloids
    alkaloid_pattern = Chem.MolFromSmarts('[#7]')  # Any nitrogen atom

    # Match using the indole skeleton pattern
    if mol.HasSubstructMatch(indole_pattern):
        # Check for nitrogen atoms indicating alkaloid characteristics
        if mol.HasSubstructMatch(alkaloid_pattern):
            return True, "Contains indole skeleton with additional alkaloid features"

        return False, "Contains indole skeleton but lacks typical alkaloid nitrogen features"
    else:
        return False, "No indole skeleton found"