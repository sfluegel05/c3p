"""
Classifies: CHEBI:28963 amino sugar
"""
from rdkit import Chem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    Amino sugars are typically sugars where a hydroxyl group is replaced by an amino group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns for detecting sugar moieties and acetamido groups
    sugar_pattern = Chem.MolFromSmarts("[C@H]1[C@@H]([O])[C@H](O)[C@H](O)[C@H](O)[C@H]1")  # Hemiacetal pyranose ring
    acetamido_pattern = Chem.MolFromSmarts("[NH]C(=O)C")  # The acetamido group structure

    # Check for sugar moiety
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar moiety found"

    # Check for acetamido group
    if not mol.HasSubstructMatch(acetamido_pattern):
        return False, "No acetamido group found"

    return True, "Contains sugar moiety with an acetamido group, indicating an amino sugar"