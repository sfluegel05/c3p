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
    
    # Look for aldehyde or ketone groups (excluding those within amide, ester, and acids)
    # Note: Redefining aldehyde/ketone SMARTS to capture more diverse cases, including aldehyde RCHO pattern
    aldehyde_ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[#6;!R]")  # Regular ketone
    aldehyde_ketone_ring_pattern = Chem.MolFromSmarts("[CX3]=[O]")  # Can be adapted to specific constraints

    # Exclude carboxylic group matches
    aldehyde_ketone_matches = mol.GetSubstructMatches(aldehyde_ketone_pattern)
    aldehyde_ketone_ring_matches = mol.GetSubstructMatches(aldehyde_ketone_ring_pattern)

    carboxylic_matches = mol.GetSubstructMatches(carboxylic_pattern)
    carboxylic_atoms = {idx for match in carboxylic_matches for idx in match}
    effective_carbonyls = [
        match for match in aldehyde_ketone_matches + aldehyde_ketone_ring_matches 
        if not set(match).intersection(carboxylic_atoms)
    ]

    # Check for at least one non-carboxylic ketone or aldehyde group
    if not effective_carbonyls:
        return False, "No distinct ketone or aldehyde group found"

    # Allow oxo fatty acids of various lengths
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:
        return False, "Not enough carbons for a fatty acid"
    
    return True, "Contains carboxylic acid group and additional oxo group(s) (aldehydic/ketonic)"