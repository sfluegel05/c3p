"""
Classifies: CHEBI:26333 prostaglandin
"""
from rdkit import Chem

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    A prostaglandin is derived from prostanoic acid with a characteristic
    cyclopentane or cyclopentene ring and specific functional groups.

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

    # SMARTS pattern for cyclopentane or cyclopentene ring
    cyclo_ring_pattern = Chem.MolFromSmarts("C1=C(C)CCC=C1 | C1=CC(=O)CCC1")
    if not mol.HasSubstructMatch(cyclo_ring_pattern):
        return (False, "No matching cyclopentane or cyclopentene ring structure found")

    # Check for terminal carboxylic acid groups
    carboxyl_group_pattern = Chem.MolFromSmarts("C(=O)O | C(=O)OC")
    if mol.GetSubstructMatches(carboxyl_group_pattern) == ():
        return (False, "No terminal carboxylic acid or ester group found")

    # Check for hydroxyl groups
    hydroxyl_group_pattern = Chem.MolFromSmarts("[C@H](O)[*]")
    if mol.GetSubstructMatch(hydroxyl_group_pattern) == ():
        return (False,"No hydroxyl groups with required stereochemistry found")

    # Check total carbon atoms, allow a slightly broader range for derivatives
    c_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if not (15 <= c_count <= 25):
        return (False, f"Unexpected carbon count: {c_count} (expected between 15 and 25)")

    return (True, "Contains key features of a prostaglandin: cyclopentane or cyclopentene ring, terminal carboxylic acid or ester, hydroxyl groups, and appropriate carbon chain length")