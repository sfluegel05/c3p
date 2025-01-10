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
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for common monosaccharide rings
    monosaccharide_patterns = [
        Chem.MolFromSmarts("[C@@H]1(CO)O[C@H](O)[C@@H](O)[C@H](O)[C@H]1O"),  # D-Glucose-like
        Chem.MolFromSmarts("[C@@H]1(CO)O[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O"),  # D-Galactose-like
        Chem.MolFromSmarts("[C@@H]1(CO)OC[C@H](O)[C@H](O)[C@H]1O"),  # D-Mannose-like
        Chem.MolFromSmarts("[C@@H]1(CO)OC[C@@H](O)[C@@H](O)[C@H](O)[C@H]1O"),  # L-fucose-like
        # Additional patterns for other saccharides can be added here
    ]
    
    # Look for rings matching saccharide substructures
    sugar_matches = sum(len(mol.GetSubstructMatches(pattern)) for pattern in monosaccharide_patterns)
    if sugar_matches != 4:
        return False, f"Only found {sugar_matches} saccharide units, need exactly 4 for a tetrasaccharide"

    # Count ether linkages as glycosidic bonds
    ether_pattern = Chem.MolFromSmarts("[C;H0]([O;H0][C;H0])")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if len(ether_matches) < 3:
        return False, "Insufficient glycosidic linkages for a tetrasaccharide"

    # Check for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 16:
        return False, "Insufficient hydroxyl groups for a tetrasaccharide"
    
    return True, "Structure matches a tetrasaccharide with four monosaccharide units linked by glycosidic bonds"