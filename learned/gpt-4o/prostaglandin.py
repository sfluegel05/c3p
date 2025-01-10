"""
Classifies: CHEBI:26333 prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    Prostaglandins are characterized by a cyclopentane ring, C20 structure, and
    presence of functional groups such as hydroxyl, carbonyl, and carboxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for cyclopentane ring (5 carbon atoms in a ring)
    cyclopentane_pattern = Chem.MolFromSmarts("C1CCCC1")
    if not mol.HasSubstructMatch(cyclopentane_pattern):
        return False, "No cyclopentane ring found"
    
    # Check for 20 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 20:
        return False, f"Incorrect number of carbon atoms: {c_count != 20}"

    # Check for hydroxyl groups - at least one
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # Check for at least one carbonyl group
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl group found"

    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OX2H]")
    carboxyl_ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    if not (mol.HasSubstructMatch(carboxyl_pattern) or mol.HasSubstructMatch(carboxyl_ester_pattern)):
        return False, "No carboxyl or ester group found"
        
    return True, "Matches characteristics of prostaglandins"