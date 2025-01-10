"""
Classifies: CHEBI:59644 oxo fatty acid
"""
from rdkit import Chem

def is_oxo_fatty_acid(smiles: str):
    """
    Determine if a molecule is an oxo fatty acid based on its SMILES string.
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

    # Look for carboxylic acid group pattern (C(=O)[O, OH])
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"

    # Look for aldehyde (R-CHO) and ketone (R2C=O) distinct from the carboxyl group
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[CH2,CH3]")
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)C")
    
    # Match patterns
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_pattern)

    # Gather carboxylic acid atoms
    carboxylic_atoms = {idx for match in carboxylic_matches for idx in match}

    # Identify legitimate oxo groups, excluding matches part of carboxylic acids
    effective_oxo_groups = [
        match for match in aldehyde_matches + ketone_matches 
        if not set(match).intersection(carboxylic_atoms)
    ]

    if not effective_oxo_groups:
        return False, "No distinct ketone or aldehyde group found"

    # Check for sufficient carbon chain length, assuming fatty acids typically have >4 carbons
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 5:
        return False, "Not enough carbons for a fatty acid"
    
    return True, "Contains carboxylic acid group and additional oxo group(s) (aldehydic/ketonic)"