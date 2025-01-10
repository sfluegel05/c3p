"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine has a piperazine-2,5-dione skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2,5-diketopiperazine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for a 2,5-diketopiperazine
    # This pattern represents a six-membered ring containing two carbonyl groups and two nitrogen atoms
    diketopiperazine_core = Chem.MolFromSmarts("O=C1NCC(=O)NC1")
    
    if mol.HasSubstructMatch(diketopiperazine_core):
        return True, "Contains a piperazine-2,5-dione skeleton"
        
    return False, "Does not contain a piperazine-2,5-dione skeleton"