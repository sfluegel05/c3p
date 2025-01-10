"""
Classifies: CHEBI:26333 prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    A prostaglandin is typically a derivative of prostanoic acid with a characteristic cyclopentane ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a 5-membered ring (cyclopentane)
    cyclopentane_pattern = Chem.MolFromSmarts("C1CCCC1")
    if not mol.HasSubstructMatch(cyclopentane_pattern):
        return False, "No cyclopentane ring found"

    # Check for carboxylic acid groups (as part of prostanoic acid backbone)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count carbon atoms to ensure structure is similar to C20 prostanoic acid
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15 or c_count > 21:  # Allowing a range around C20
        return False, f"Unexpected carbon count: {c_count} (expected ~20)"

    # Check for at least one hydroxyl group (-OH), which is common in prostaglandins
    hydroxyl_pattern = Chem.MolFromSmarts("CO")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # The minimum required features are present, classify as prostaglandin
    return True, "Contains key features of a prostaglandin: cyclopentane ring, carboxylic acid, and hydroxyl groups"