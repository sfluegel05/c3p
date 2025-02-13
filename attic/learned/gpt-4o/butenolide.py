"""
Classifies: CHEBI:50523 butenolide
"""
from rdkit import Chem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide is defined as a gamma-lactone with a 2-furanone skeleton and its substituted derivatives.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general SMARTS pattern for the butenolide core (gamma-lactone)
    # allowing for flexibility in substitutions at various positions
    butenolide_pattern = Chem.MolFromSmarts("O=C1O[C;!H1]C=C1")  # General gamma-lactone
    
    # Additional patterns to capture variations
    extra_patterns = [
        Chem.MolFromSmarts("O=C1O[C;!H1]=CC1"),     # Double bond within ring
        Chem.MolFromSmarts("O=C1O[C;!H1]C=C1[!c]"), # Allow exocyclic unsaturation
        Chem.MolFromSmarts("O=C1C=CC(=O)O1"),       # Account for oxo and hydroxyl groups
        Chem.MolFromSmarts("O=C1C=COC1"),           # Flexibility in double bond location
    ]
    
    # Check for the base substructure match
    if mol.HasSubstructMatch(butenolide_pattern):
        return True, "Contains 2-furanone skeleton typical of a butenolide"

    # Check for extra pattern matches
    for pattern in extra_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a variant of the 2-furanone/gamma-lactone skeleton"

    return False, "Does not contain a 2-furanone skeleton or any recognized substituted form"