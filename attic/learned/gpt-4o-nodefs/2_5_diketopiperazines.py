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
    
    # Define more comprehensive SMARTS patterns to match 2,5-diketopiperazine variability
    # Pattern accounts for possible substitutions and flexibility around the cyclic core
    diketopiperazine_patterns = [
        Chem.MolFromSmarts("N1C(=O)CNC(=O)C1"), # Simple diketopiperazine
        Chem.MolFromSmarts("N1C(=O)[C@H,R][NR][C@H,R]C(=O)1"), # Include stereochemistry
        Chem.MolFromSmarts("O=C1N[C@H,R]C(=O)N[C@H,R]1"), # Flexible connections
    ]
    
    for pattern in diketopiperazine_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains 2,5-diketopiperazine motif"
    
    return False, "Does not contain 2,5-diketopiperazine motif"