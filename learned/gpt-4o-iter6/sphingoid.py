"""
Classifies: CHEBI:35785 sphingoid
"""
from rdkit import Chem

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    Sphingoids are primarily characterized by a long aliphatic chain with one or more
    hydroxyl groups and an amino group, with possible unsaturation.

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
    
    # Refined pattern for detecting long aliphatic chain
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC[CH3;CH2]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Does not have a long enough aliphatic carbon chain"

    # Refined amino group pattern - includes primary and secondary amines and associated functionality
    amino_pattern = Chem.MolFromSmarts("[NX3][CX4H2]")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "Missing amino group"
    
    # Refine hydroxyl pattern checking especially for vicinal arrangement (important in sphingoids)
    hydroxyl_pattern = Chem.MolFromSmarts("[CX3,CX4][OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing hydroxyl group"
    
    # Check for optional unsaturations (double bonds), while not required validation as sphingoids can be fully saturated
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    
    # Detection logic to classify
    if mol.HasSubstructMatch(double_bond_pattern):
        return True, "Contains long aliphatic chain with hydroxyl and amino groups, including possible unsaturation; could be a sphingoid."
    else:
        return True, "Contains long aliphatic chain with hydroxyl and amino groups; even if fully saturated, it could be a sphingoid."