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

    # Identify saccharide units - presence of multiple cyclic ethers with stereochemistry and hydroxyl groups
    saccharide_pattern = Chem.MolFromSmarts("[C@H]1([O,R])[C@H]([O,R])[C@H]([O,R])C([O,R])[C@H]1O")
    saccharide_matches = mol.GetSubstructMatches(saccharide_pattern)

    if len(saccharide_matches) < 4:
        return False, f"Less than 4 monosaccharide units found, found {len(saccharide_matches)}"

    # Verify glycosidic linkages - Oxygen bridges between rings
    linkage_pattern = Chem.MolFromSmarts("O[C@H]1[C@H]([O,R])[C@H]([O,R])C([O,R])[C@H]1O")
    linkage_matches = mol.GetSubstructMatches(linkage_pattern)

    if len(linkage_matches) < 3:
        return False, f"Insufficient glycosidic linkages, found {len(linkage_matches)}"

    return True, "Molecule is a tetrasaccharide with appropriate saccharide units and glycosidic linkages"