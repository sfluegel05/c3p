"""
Classifies: CHEBI:26333 prostaglandin
"""
from rdkit import Chem

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    Prostaglandins are typically derived from the parent C20 prostanoic acid,
    featuring a cyclopentane ring with various functional modifications.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a basic cyclopentane ring pattern with flexibility for substitutions
    cyclopentane_pattern = Chem.MolFromSmarts("C1CCC(C)C1")
    if cyclopentane_pattern is None:  # Ensure the pattern is valid
        return None, "Internal Error: cyclopentane pattern could not be created"
    
    if not mol.HasSubstructMatch(cyclopentane_pattern):
        return False, "No cyclopentane ring structure found"

    # Check carbon count; should be around 20 for prostaglandins
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (18 <= c_count <= 22):
        return False, f"Carbon count out of range (expected around 20, got {c_count})"

    # Check for at least one hydroxyl group - an essential feature
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if hydroxyl_pattern is None:  # Ensure pattern validity
        return None, "Internal Error: hydroxyl pattern could not be created"
    
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # Check for at least one carbonyl group - another essential feature
    carbonyl_pattern = Chem.MolFromSmarts("C=O")
    if carbonyl_pattern is None:
        return None, "Internal Error: carbonyl pattern could not be created"
    
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl group found"

    # Check for a carboxylic acid or ester group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OX2H]")
    if carboxyl_pattern is None:
        return None, "Internal Error: carboxyl pattern could not be created"

    carboxyl_ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    if carboxyl_ester_pattern is None:
        return None, "Internal Error: carboxyl/ester pattern could not be created"
    
    if not (mol.HasSubstructMatch(carboxyl_pattern) or mol.HasSubstructMatch(carboxyl_ester_pattern)):
        return False, "No carboxylic acid or ester group found"

    return True, "Contains features characteristic of prostaglandins"