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

    # Look for a 5-membered ring with possible unsaturation (e.g., cyclopentenone)
    cyclopentane_pattern = Chem.MolFromSmarts("C1=CC=CC1")
    if not mol.HasSubstructMatch(cyclopentane_pattern):
        return False, "No cyclopentene/cyclopentane ring found"

    # Check for carboxylic acid or ester groups
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid or ester group found"

    # Verify presence of hydroxyl group (-OH) or other oxygen-containing groups suggesting biological activity
    hydroxyl_or_ether_pattern = Chem.MolFromSmarts("CO")
    if not mol.HasSubstructMatch(hydroxyl_or_ether_pattern):
        return False, "No hydroxyl or ether group found"

    # Slightly relax carbon count to allow for typical modifications
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15 or c_count > 25:  # Allow more flexibility
        return False, f"Unexpected carbon count: {c_count} (expected ~20)"

    # The minimum required features are present, classify as prostaglandin
    return True, "Contains key features of a prostaglandin: a flexible ring structure, carboxylic acid/ester, and hydroxyl/ether groups"