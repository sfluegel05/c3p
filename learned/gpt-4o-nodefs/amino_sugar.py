"""
Classifies: CHEBI:28963 amino sugar
"""
from rdkit import Chem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    Amino sugars are typically sugars where a hydroxyl group is replaced by an amino group,
    often with an acetamido modification.

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

    # Define patterns for detecting sugar moieties and amino group modifications
    # A common pattern in amino sugars would be a cyclic structure with several hydroxyls
    sugar_pattern = Chem.MolFromSmarts("[C@H]1([O])[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@@H]1")  # Pyranose-like ring
    acetamido_pattern = Chem.MolFromSmarts("N[C@H]C(=O)C")  # Common pattern for N-acetyl group

    # Check for sugar moiety
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar moiety found"

    # Check for acetamido group (or amino group in general)
    if not mol.HasSubstructMatch(acetamido_pattern):
        # If the acetamido is not found, try a basic amino check
        if mol.HasSubstructMatch(Chem.MolFromSmarts("N")):
            return True, "Contains sugar moiety with an amino group, indicating an amino sugar"
        else:
            return False, "No acetamido or amino group modification found"

    return True, "Contains sugar moiety with an acetamido group, indicating an amino sugar"