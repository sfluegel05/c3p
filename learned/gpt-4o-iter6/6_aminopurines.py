"""
Classifies: CHEBI:20706 6-aminopurines
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule belongs to the 6-aminopurines class based on its SMILES string.
    6-aminopurine (adenine) is characterized by a specific purine structure with an amino group at the 6th position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains 6-aminopurine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define adenine substructure patterns with improved coverage
    adenine_patterns = [
        Chem.MolFromSmarts("n1cnc2c(nc[nH]c2)n1"),  # Basic adenine pattern
        Chem.MolFromSmarts("c1[nH]c2c(ncnc2)n1"),  # Another common structure
        Chem.MolFromSmarts("n1cnc2c(nc[nH]c2)n1"), # Variant with different hydrogen positions
        Chem.MolFromSmarts("n1c[nH]c2cncnc12"),    # Recognition of other tautomeric variations
        Chem.MolFromSmarts("n1cnc2[nH]cnc2c1"),    # Pattern for common protonation states
    ]
    
    # Search for any adenine substructure in the molecule
    for pattern in adenine_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains 6-aminopurine (adenine) as part of its structure"
    
    return False, "Does not contain 6-aminopurine (adenine) structure"