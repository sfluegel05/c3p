"""
Classifies: CHEBI:28963 amino sugar
"""
from rdkit import Chem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    Amino sugars are typically sugars where a hydroxyl group is replaced by an amino group,
    often with an N-acetamido modification.

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

    # Define patterns for sugar rings and amino group modifications
    # Assume pyranose ring common for sugars, with at least one amino or acetamido replacement
    sugar_pattern = Chem.MolFromSmarts("[C@H]1(O)[C@@H]([O])[C@H](O)[C@@H](O)[C@@H]([O])[C@@H]1")  # General pyranose-like
    amino_pattern = Chem.MolFromSmarts("[C@H]1(O)[C@@H](N)[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1")  # Pyranose with amino
    acetamido_pattern = Chem.MolFromSmarts("N[C@H](C)C(=O)")  # Common pattern for N-acetyl group

    # Check for sugar moiety
    has_sugar = mol.HasSubstructMatch(sugar_pattern) or mol.HasSubstructMatch(amino_pattern)
    if not has_sugar:
        return False, "No sugar moiety found"

    # Check for either amino group or acetylation
    has_amino_modification = mol.HasSubstructMatch(amino_pattern) or mol.HasSubstructMatch(acetamido_pattern)
    if not has_amino_modification:
        return False, "No amino or acetamido group modification found"

    return True, "Contains sugar moiety with appropriate amino group modification, indicating an amino sugar"