"""
Classifies: CHEBI:50126 tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Define SMARTS patterns for common monosaccharide units
    common_sugar_patterns = [
        Chem.MolFromSmarts("C1(O)CO[C@@H](O)[C@H](O)C1"),  # For example, D-glucose
        # Add more specific patterns here for other sugars if needed
    ]
    
    # Look for rings matching saccharide substructures
    sugar_matches = sum(mol.HasSubstructMatch(pattern) for pattern in common_sugar_patterns)
    if sugar_matches < 4:
        return False, f"Only found {sugar_matches} saccharide units, need exactly 4 for a tetrasaccharide"

    # Count the ether linkages as glycosidic bonds
    ether_pattern = Chem.MolFromSmarts("C-O-C")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if len(ether_matches) < 3:
        return False, "Insufficient glycosidic linkages for a tetrasaccharide"

    # Check for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 8:
        return False, "Insufficient hydroxyl groups for a tetrasaccharide"
    
    return True, "Structure matches a tetrasaccharide with four monosaccharide units linked by glycosidic bonds"