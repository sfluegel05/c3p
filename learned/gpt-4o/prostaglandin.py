"""
Classifies: CHEBI:26333 prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    Prostaglandins are characterized by a core cyclopentane or similar ring,
    and presence of functional groups such as hydroxyl, carbonyl, and carboxyl groups.
    Variations in carbon count are allowed based on known derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible pattern for a cyclopentane or related ring (e.g., including double bonds)
    ring_pattern = Chem.MolFromSmarts("C1=CC=CC1 | C1CCC(C=O)C1")
    if not mol.HasSubstructMatch(ring_pattern):
        return False, "No cyclopentane-like ring found"

    # Accept a range in carbon count (e.g., 18-22 to accommodate variations)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (18 <= c_count <= 22):
        return False, f"Carbon count out of range (expected 18-22, got {c_count})"
    
    # Check for hydroxyl groups - at least one
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # Check for at least one carbonyl group
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl group found"

    # Check for carboxylic acid or ester group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OX2H]")
    carboxyl_ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    if not (mol.HasSubstructMatch(carboxyl_pattern) or mol.HasSubstructMatch(carboxyl_ester_pattern)):
        return False, "No carboxyl or ester group found"
        
    return True, "Matches characteristics of prostaglandins"