"""
Classifies: CHEBI:20706 6-aminopurines
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule contains the 6-aminopurine structure
    (adenine) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains the 6-aminopurine structure,
              False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for 6-aminopurine structure: adenine with an amine at position 6
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(N)nc[nH]c2n1")  # Correct adenine pattern with amine at position 6

    if adenine_pattern is None:
        return False, "Failed to create 6-aminopurine pattern"
    
    # Check for the adenine substructure with the amine group
    if mol.HasSubstructMatch(adenine_pattern):
        return True, "Contains 6-aminopurine (adenine) structure"
    
    return False, "Does not contain 6-aminopurine (adenine) structure"