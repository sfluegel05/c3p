"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine is characterized by a six-membered ring with two nitrogen atoms
    and ketone groups on the carbons adjacent to these nitrogens.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2,5-diketopiperazine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved SMARTS pattern for 2,5-diketopiperazines:
    # A six-membered ring (with any carbon substitution) containing C(=O)N and NC(=O)
    # This pattern ensures only the core structure is specified with flexibility in substitutions
    diketopiperazine_pattern = Chem.MolFromSmarts("N1C(=O)CNC(=O)C1")

    if mol.HasSubstructMatch(diketopiperazine_pattern):
        return True, "Contains a piperazine-2,5-dione skeleton"
    
    return False, "Does not contain a piperazine-2,5-dione skeleton"