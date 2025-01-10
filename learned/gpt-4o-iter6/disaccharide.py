"""
Classifies: CHEBI:36233 disaccharide
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is two monosaccharides joined by a glycosidic bond.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns for 5- and 6-membered sugar rings (furanose and pyranose)
    sugar_patterns = [
        Chem.MolFromSmarts("O1[C@@H](CO)[C@H](O)[C@@H](O)[C@H]1"),  # 5-Membered ring
        Chem.MolFromSmarts("O1[C@H](O)[C@H](O[C@@H](O)C)[C@@H](O)C1"),  # 6-Membered ring (example)
        Chem.MolFromSmarts("O1[C@H](O[C@@H](O)[C@H](O)[C@@]1)C"),  # Generic 6-membered, reverse roles
    ]

    # Detect sugar rings
    sugar_ring_count = 0
    for pattern in sugar_patterns:
        matches = mol.GetSubstructMatches(pattern)
        sugar_ring_count += len(matches)

    if sugar_ring_count < 2:
        return False, f"Expected 2 sugar rings, found {sugar_ring_count}"

    # Patterns for glycosidic linkage (anomeric carbon connected through oxygen)
    glycosidic_patterns = [
        Chem.MolFromSmarts("[C@H]1(O[C@H]2[*:1])[*:2]1C2"),
        Chem.MolFromSmarts("[C@H]1(O[C@@H]2[*:1])[*:2]1C2"),
        Chem.MolFromSmarts("[C@@H]1(O[C@H]2[*:1])[*:2]1C2"),
        Chem.MolFromSmarts("[C@@H]1(O[C@@H]2[*:1])[*:2]1C2"),
    ]

    # Check for at least one valid glycosidic linkage
    glyco_bond_count = 0
    for pattern in glycosidic_patterns:
        matches = mol.GetSubstructMatches(pattern)
        glyco_bond_count += len(matches)

    if glyco_bond_count < 1:
        return False, "No appropriate glycosidic link found to connect sugars"

    return True, "Contains two sugar rings joined by a glycosidic bond"