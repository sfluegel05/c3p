"""
Classifies: CHEBI:50126 tetrasaccharide
"""
from rdkit import Chem

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide comprises four monosaccharide units linked by glycosidic bonds.

    Args:
        smiles (str): SMILES representation of the molecule

    Returns:
        bool: True if the molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define extended SMARTS patterns for common monosaccharide rings, considering stereochemistry and oxygen connectivity in sugars
    monosaccharide_patterns = [
        Chem.MolFromSmarts("[C@H]1(O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)O1"),  # D-Glucose-like
        Chem.MolFromSmarts("[C@H]1(O)[C@H](O)[C@H](O)[C@H](O)[C@@H](O)O1"),  # D-Galactose-like
        Chem.MolFromSmarts("[C@H]1(O)[C@H](O)[C@@H](O)[C@H](O)O1"),  # D-Mannose-like
        Chem.MolFromSmarts("[C@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)O1"),  # L-fucose/furanose-like
        # More patterns can be added for completeness
    ]
    
    sugar_matches = sum(len(mol.GetSubstructMatches(pattern)) for pattern in monosaccharide_patterns)
    if sugar_matches != 4:
        return False, f"Only found {sugar_matches} saccharide units, need exactly 4 for a tetrasaccharide"

    # Count ether linkages, also considering the O-C-O pattern which is characteristic of glycosidic bonds
    ether_pattern = Chem.MolFromSmarts("[C;H0][OX2][C;H0]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if len(ether_matches) < 3:
        return False, f"Insufficient glycosidic linkages for a tetrasaccharide, found {len(ether_matches)}"

    # Check for an adequate number of hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 8 * sugar_matches // 2:
        return False, f"Insufficient hydroxyl groups for a tetrasaccharide, found {len(hydroxyl_matches)}"
    
    return True, "Structure matches a tetrasaccharide with four monosaccharide units linked by glycosidic bonds"