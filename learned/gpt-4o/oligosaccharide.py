"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is an oligosaccharide, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for common hexose sugar rings (like glucose, galactose, mannose)
    sugar_patterns = [
        Chem.MolFromSmarts("[C@H]1(O)[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"), # Hexopyranose pattern
        Chem.MolFromSmarts("[C@@H]1(O)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O")  # Mirror stereochemistry
    ]
    
    # Define SMARTS pattern for glycosidic linkages (ether links between sugars)
    glycosidic_bond_pattern = Chem.MolFromSmarts("O[C@]")  # Example ether linkage pattern

    # Check for at least two sugar units
    sugar_count = 0
    for pattern in sugar_patterns:
        sugar_count += len(mol.GetSubstructMatches(pattern))
    
    if sugar_count < 2:
        return False, "Insufficient monosaccharide units"

    # Check for glycosidic linkages connecting sugars
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic linkages found"

    return True, "Contains sufficient monosaccharide units and glycosidic linkages indicative of an oligosaccharide"

# This method improves on exact matches for oligosaccharides by using common sugar ring structures and checking for glycosidic linkages.