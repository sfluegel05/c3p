"""
Classifies: CHEBI:26333 prostaglandin
"""
from rdkit import Chem

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    A prostaglandin is a derivative of prostanoic acid with characteristic cyclopentane
    or cyclopentene core and specific functional groups indicating biological activity.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")

    # Improved SMARTS patterns for cyclopentane and cyclopentene with specific stereo configurations
    cyclopentane_pattern = Chem.MolFromSmarts("[C@H]1[C@H]([C@H]([C@H]([C@H]1)*)*)*")
    cyclopentene_pattern = Chem.MolFromSmarts("[C@H]1[C@H]([C@H]([C@H]([C@@H]1*)*)*)*")
    
    # Detection of cyclopentane/cyclopentene with stereo configurations
    has_cyclo_ring = mol.HasSubstructMatch(cyclopentane_pattern) or mol.HasSubstructMatch(cyclopentene_pattern)
    if not has_cyclo_ring:
        return (False, "No matching cyclopentane or cyclopentene ring structure found")

    # Refined pattern for carboxylic acid and ester groups
    carboxyl_group_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_group_pattern = Chem.MolFromSmarts("C(=O)OC")
    
    if not (mol.HasSubstructMatch(carboxyl_group_pattern) or mol.HasSubstructMatch(ester_group_pattern)):
        return (False, "No carboxylic acid or ester group found")

    # Look for hydroxyl groups with specified stereo-configuration, typical for prostaglandins
    hydroxyl_group_pattern = Chem.MolFromSmarts("[C@H](O)[C@H]") 
    ether_group_pattern = Chem.MolFromSmarts("C-O-C")
    
    has_hydroxyl_or_ether = mol.HasSubstructMatch(hydroxyl_group_pattern) or mol.HasSubstructMatch(ether_group_pattern)
    if not has_hydroxyl_or_ether:
        return (False, "No hydroxyl or ether group found")

    # Check for a suitable chain length typical for prostaglandins
    c_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if not (15 <= c_count <= 30):
        return (False, f"Unexpected carbon count: {c_count} (expected between 15 and 30)")

    # If all checks pass, classify it as a prostaglandin
    return (True, "Contains key features of a prostaglandin: characteristic ring structure, carboxyl/ester groups, and hydroxyl groups")