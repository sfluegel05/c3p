"""
Classifies: CHEBI:20706 6-aminopurines
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule is a 6-aminopurine (adenine derivative) based on its SMILES string.
    A 6-aminopurine (adenine) has a purine base with an amino group at the 6-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 6-aminopurine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the 6-aminopurine (adenine) structure
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)N")
    if mol.HasSubstructMatch(adenine_pattern):
        return True, "Contains adenine structure"

    return False, "No adenine structure found"