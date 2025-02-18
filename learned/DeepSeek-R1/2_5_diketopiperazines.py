"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: 2,5-diketopiperazines (CHEBI:53758)
"""
from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine has a piperazine-2,5-dione skeleton (six-membered ring with two ketone groups at positions 2 and 5).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2,5-diketopiperazine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for the piperazine-2,5-dione core
    dkp_pattern = Chem.MolFromSmarts("[O]=C1NCC(=O)NC1")
    
    # Check if the pattern is present in the molecule
    if mol.HasSubstructMatch(dkp_pattern):
        return True, "Contains piperazine-2,5-dione skeleton"
    else:
        return False, "No piperazine-2,5-dione skeleton found"