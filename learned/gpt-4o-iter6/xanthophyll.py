"""
Classifies: CHEBI:27325 xanthophyll
"""
from rdkit import Chem

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    A xanthophyll is a carotenoid that is oxygenated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # General pattern representing alternation of double and single bonds (conjugation common in carotenoids)
    conjugation_pattern = Chem.MolFromSmarts("C=C(-,=C)*=C")
    if conjugation_pattern and not mol.HasSubstructMatch(conjugation_pattern):
        return False, "No proper alternating double-bonds pattern typically found in carotenoids"
    
    # Check for the presence of oxygen atoms in various functional groups
    oxygen_pattern = Chem.MolFromSmarts("[OX1,OX2,OX3]")
    if oxygen_pattern and not mol.HasSubstructMatch(oxygen_pattern):
        return False, "Carotenoid derivative must be oxygenated (No oxygen detected)"
    
    # Check for possible multiple types of oxygen functionalization present
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")  # Detect -OH groups
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")  # Detect =O groups
    ether_pattern = Chem.MolFromSmarts("[CX3][OX2R]")  # Detect -O- groups

    if not (mol.HasSubstructMatch(hydroxyl_pattern) or mol.HasSubstructMatch(carbonyl_pattern) or mol.HasSubstructMatch(ether_pattern)):
        return False, "At least one oxygen functional group required (such as hydroxyl, carbonyl, or ether not found)"

    # Additional pattern checks and logic improvements as needed

    return True, "Contains features of an oxygenated carotenoid (xanthophyll)"