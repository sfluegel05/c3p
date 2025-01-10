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

    # Look for ether groups: indicative of glycosidic linkages
    ether_pattern = Chem.MolFromSmarts("C-O-C")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if len(ether_matches) < 3:  # Expect at least 3 ether linkages for 4 monosaccharides
        return False, "Insufficient glycosidic linkages for tetrasaccharide"

    # Look for hydroxyl groups: common in sugar units
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 8:  # Estimate at least 2 OH per monosaccharide
        return False, "Insufficient hydroxyl groups for a tetrasaccharide"

    # Check for ring structures typical in saccharides
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 4:
        return False, "Insufficient ring structures for four saccharide units"

    # Check if the molecule consists of 4 saccharide rings using specific SMARTS pattern for sugars
    # This SMARTS pattern might match common monosaccharide ring structures
    sugar_ring_pattern = Chem.MolFromSmarts("C1(O)CO[C@@H](O)[C@H](O)C1")
    sugar_ring_matches = mol.GetSubstructMatches(sugar_ring_pattern)
    if len(sugar_ring_matches) != 4:
        return False, f"Incorrect number of sugar rings ({len(sugar_ring_matches)}) for tetrasaccharide"
    
    return True, "Structure matches a tetrasaccharide with four monosaccharide units linked by glycosidic bonds"