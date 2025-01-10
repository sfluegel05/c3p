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

    # Define a broader SMARTS pattern for the butenolide core
    # This accounts for the base 2-furanone and allows various substitutions
    butenolide_pattern = Chem.MolFromSmarts("O=C1OC=CC1")  # Base furanone pattern
    extra_patterns = [
        Chem.MolFromSmarts("O=C1C=COC1"),  # Account for different unsaturation
        Chem.MolFromSmarts("O=C1OC=CC1=O"),  # Account for further oxidation
        Chem.MolFromSmarts("O=C1C=CC(=O)O1")  # Allow for broader keto-lactone variability
    ]
    
    # Check for the base substructure match
    if mol.HasSubstructMatch(butenolide_pattern):
        return True, "Contains 2-furanone skeleton typical of a butenolide"

    # Check for extra pattern matches
    for pattern in extra_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a variant of the 2-furanone/gamma-lactone skeleton"

    return False, "Does not contain a 2-furanone skeleton or any recognized substituted form"