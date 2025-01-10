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
    diketopiperazine_pattern = Chem.MolFromSmarts("C1=NC(=O)CC(=O)N1")
    if mol.HasSubstructMatch(diketopiperazine_pattern):
        return True, "Contains 2,5-diketopiperazine motif"
    
    return False, "Does not contain 2,5-diketopiperazine motif"