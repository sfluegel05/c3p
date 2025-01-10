"""
Classifies: CHEBI:51963 hopanoid
"""
from rdkit import Chem

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    Hopanoids are a class of pentacyclic triterpenoids with hopane backbones.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # A more generalized pattern to capture hopanoid backbones, including functional variations
    hopanoid_patterns = [
        Chem.MolFromSmarts("C1CCC2CC3CC4C5CC(C4)C5CC3C2C1"),  # Basic hopane pattern
        Chem.MolFromSmarts("C1CCC2CC3CC(C4)C4CC3C2C1"),        # Allow for modified hopanes
        Chem.MolFromSmarts("C1CCC2C(C3CC(C4)C4CCC23)CCC1"),     # Other cyclization variations
    ]

    for hopanoid_pattern in hopanoid_patterns:
        if mol.HasSubstructMatch(hopanoid_pattern):
            return True, "Contains hopanoid-like backbone"

    return False, "No hopanoid backbone found"

# This updated pattern attempts to broadly detect the triterpenoid backbone characteristic of hopanoids.