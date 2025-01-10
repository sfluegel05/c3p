"""
Classifies: CHEBI:50126 tetrasaccharide
"""
from rdkit import Chem

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide is a carbohydrate composed of four saccharide units linked by glycosidic bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define general sugar ring pattern
    cyclic_sugar_pattern = Chem.MolFromSmarts("OC1CO[C@H]([C@H](O)[C@@H](O)[C@H]1O)")  # generic sugar ring
    # Define glycosidic bond pattern
    glycosidic_linkage_pattern = Chem.MolFromSmarts("-O-")

    # Initialize counts
    saccharide_count = 0
    linkage_count = 0

    # Check for saccharide units using a generic sugar pattern
    saccharide_matches = mol.GetSubstructMatches(cyclic_sugar_pattern)
    saccharide_count = len(saccharide_matches)
    
    # Check for glycosidic linkages, more generic linkage pattern
    linkage_matches = mol.GetSubstructMatches(glycosidic_linkage_pattern)
    linkage_count = len(linkage_matches)

    # Assess if it's a tetrasaccharide (assuming linear unbranched)
    if saccharide_count >= 4 and linkage_count >= 3:
        return True, "Molecule is a tetrasaccharide with appropriate saccharide units and glycosidic linkages"

    # Provide reason for failure if not a tetrasaccharide
    if saccharide_count < 4:
        return False, f"Less than 4 monosaccharide units found, found {saccharide_count}"
    if linkage_count < 3:
        return False, f"Insufficient glycosidic linkages, found {linkage_count}"

    return False, "Unclassified - unexpected structure"