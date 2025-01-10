"""
Classifies: CHEBI:26333 prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Define SMARTS pattern for cyclopentane or cyclopentene ring
    cyclopentane_pattern = Chem.MolFromSmarts("[CX4]1[CX4][CX4][CX4][CX4]1")
    cyclopentene_pattern = Chem.MolFromSmarts("[C,c]1([C,c])[C,c][C,c][C,c]1")

    # Determine if the molecule has a cyclopentane or cyclopentene ring
    has_cyclo_ring = mol.HasSubstructMatch(cyclopentane_pattern) or mol.HasSubstructMatch(cyclopentene_pattern)
    if not has_cyclo_ring:
        return (False, "No cyclopentane or cyclopentene ring found")

    # Look for carboxylic acid or ester functional group
    carboxylic_acid_or_ester_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(carboxylic_acid_or_ester_pattern):
        return (False, "No carboxylic acid or ester group found")

    # Ensure presence of hydroxyl group (-OH) or ether linkage
    hydroxyl_or_ether_pattern = Chem.MolFromSmarts("[OX2H][CX4]")  # Matches C-O-H or C-O-C
    if not mol.HasSubstructMatch(hydroxyl_or_ether_pattern):
        return (False, "No hydroxyl or ether group found")

    # Verify carbon count
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15 or c_count > 30:  # Allow more flexibility, prostaglandins have around 20 carbons
        return (False, f"Unexpected carbon count: {c_count} (expected ~20)")

    # The minimum required features are present, classify as prostaglandin
    return (True, "Contains key features of a prostaglandin: essential ring structure, carboxylic/ester group, and hydroxyl/ether groups")