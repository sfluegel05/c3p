"""
Classifies: CHEBI:20706 6-aminopurines
"""
"""
Classifies: CHEBI:26334 6-aminopurines
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule is a 6-aminopurine based on its SMILES string.
    A 6-aminopurine is any compound having 6-aminopurine (adenine) as part of its structure.

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

    # Define the adenine substructure pattern
    adenine_pattern = Chem.MolFromSmarts("N1C=NC2=C1N=CN=C2N")
    
    # Check if the molecule contains the adenine substructure
    if mol.HasSubstructMatch(adenine_pattern):
        return True, "Contains 6-aminopurine (adenine) substructure"
    else:
        return False, "Does not contain 6-aminopurine (adenine) substructure"