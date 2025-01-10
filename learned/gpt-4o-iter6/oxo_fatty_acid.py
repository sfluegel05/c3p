"""
Classifies: CHEBI:59644 oxo fatty acid
"""
from rdkit import Chem

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid is defined as any fatty acid containing at least one aldehydic or ketonic group
    in addition to the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group pattern (C(=O)O)
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for distinct aldehyde (R-CHO) and ketone (R2C=O) groups
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[CX4,CH2,CH3]")  # Aldehyde R-CHO
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[CX4]")  # Ketone R2C=O
    
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)

    # Ensure those matches are not part of carboxylic acids
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_pattern)
    carboxylic_atoms = {idx for match in carboxylic_matches for idx in match}
    
    # Identify effective oxo groups, excluding carboxyl carbon and oxygen
    effective_oxo_groups = [
        match for match in aldehyde_matches + ketone_matches 
        if not set(match).intersection(carboxylic_atoms)
    ]

    if not effective_oxo_groups:
        return False, "No distinct ketone or aldehyde group found"

    # Check the length of the carbon chain to confirm fatty acid nature
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:
        return False, "Not enough carbons for a fatty acid"
    
    return True, "Contains carboxylic acid group and additional oxo group(s) (aldehydic/ketonic)"