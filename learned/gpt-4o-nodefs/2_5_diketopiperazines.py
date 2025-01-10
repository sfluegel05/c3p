"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2,5-diketopiperazine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for a 2,5-diketopiperazine
    # This pattern is more flexible and accounts for typical diketo structures
    diketopiperazine_pattern = Chem.MolFromSmarts("O=C1NCC(=O)N1")
    
    if mol.HasSubstructMatch(diketopiperazine_pattern):
        # Further checks can be added here if necessary to refine classification
        return True, "Contains 2,5-diketopiperazine motif"
    
    # Additional exclusions if necessary (e.g., to exclude barbiturate-like structures)
    
    return False, "Does not contain 2,5-diketopiperazine motif"