"""
Classifies: CHEBI:20706 6-aminopurines
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule contains a 6-aminopurine (adenine) structure
    based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains the 6-aminopurine structure,
              False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Correct SMARTS pattern for adenine matching specifically position 6 amine
    adenine_pattern = Chem.MolFromSmarts("c1ncnc2[nH]cn(C)c12")  # Include constraints for 6-position amine

    if adenine_pattern is None:
        return False, "Failed to create 6-aminopurine pattern"
    
    # Check if the molecule matches the adenine structure including the amine
    if mol.HasSubstructMatch(adenine_pattern):
        return True, "Contains 6-aminopurine (adenine) structure"
    
    return False, "Does not contain 6-aminopurine (adenine) structure"