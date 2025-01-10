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

    # Pattern for cyclopentane or similar ring (flexible, with double bonds)
    cyclopentane_pattern = Chem.MolFromSmarts("C1CCCC1 | C1CCC(C=O)C1 | C1=CC=CC1")
    if not mol.HasSubstructMatch(cyclopentane_pattern):
        return False, "No cyclopentane or similar ring structure found"

    # Check for carbon count around typical range for prostaglandins (18-22 as a reasonable assumption)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (18 <= c_count <= 22):
        return False, f"Carbon count out of range (expected 18-22, got {c_count})"

    # Check for at least one hydroxyl group - an essential feature
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # Check for at least one carbonyl group - another essential feature
    carbonyl_pattern = Chem.MolFromSmarts("C=O")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl group found"

    # Check for a carboxylic acid or ester group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OX2H]")
    carboxyl_ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    if not (mol.HasSubstructMatch(carboxyl_pattern) or mol.HasSubstructMatch(carboxyl_ester_pattern)):
        return False, "No carboxylic acid or ester group found"

    return True, "Contains features characteristic of prostaglandins"