"""
Classifies: CHEBI:35785 sphingoid
"""
from rdkit import Chem

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    Sphingoids are primarily characterized by a long aliphatic chain with one or more
    hydroxyl groups and an amino group, with possible unsaturations.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a sphingoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for long aliphatic chain
    carbon_chain_pattern = Chem.MolFromSmarts("C(CCCCCCCCCC)CCCC")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Does not have a long aliphatic carbon chain"

    # Check for hydroxyl group(s)
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H1]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, "Missing hydroxyl group(s)"

    # Check for an amino group
    amino_pattern = Chem.MolFromSmarts("[NX3][CX4]")
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    if len(amino_matches) < 1:
        return False, "Missing amino group"
    
    # Check for optional unsaturations (double bonds)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    
    # A sphingoid may or may not contain a double bond
    if not mol.HasSubstructMatch(double_bond_pattern):
        return True, "Contains long aliphatic chain with hydroxyl and amino groups; lacks unsaturation, but still could be a sphingoid."
    
    return True, "Contains long aliphatic chain with hydroxyl and amino groups, including possible unsaturation; could be a sphingoid"