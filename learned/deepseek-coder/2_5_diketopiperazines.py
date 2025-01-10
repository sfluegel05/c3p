"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: CHEBI:47778 2,5-diketopiperazine
"""
from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine has a piperazine-2,5-dione skeleton.

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

    # Define the piperazine-2,5-dione skeleton pattern
    # The pattern matches a six-membered ring with two nitrogen atoms at positions 1 and 4,
    # and two carbonyl groups at positions 2 and 5.
    piperazine_dione_pattern = Chem.MolFromSmarts("[O]=C1NCC(=O)NC1")
    
    # Check if the molecule contains the piperazine-2,5-dione skeleton
    if not mol.HasSubstructMatch(piperazine_dione_pattern):
        return False, "No piperazine-2,5-dione skeleton found"

    return True, "Contains a piperazine-2,5-dione skeleton"