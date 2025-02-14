"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide is an oligosaccharide comprising four monomeric monosaccharide units.

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

    # Define a SMARTS pattern for monosaccharide rings (both 5 and 6 membered rings)
    # The carbon with two oxygens is explicit, as it is part of the cycle of monosaccharides
    # Include all possible ring configurations with explicit CH, CH2 and oxygen atoms
    monosaccharide_pattern_6 = Chem.MolFromSmarts("[CX4H]([OX2])[CX4H][CX4H][CX4H][CX4H]([OX2])1")
    monosaccharide_pattern_5 = Chem.MolFromSmarts("[CX4H]([OX2])[CX4H][CX4H][CX4H]([OX2])1")
    if monosaccharide_pattern_6 is None or monosaccharide_pattern_5 is None:
        return None, "Error in SMARTS pattern"

    # Find all matching monosaccharide units
    monosaccharide_matches_6 = mol.GetSubstructMatches(monosaccharide_pattern_6)
    monosaccharide_matches_5 = mol.GetSubstructMatches(monosaccharide_pattern_5)
    monosaccharide_matches = monosaccharide_matches_6 + monosaccharide_matches_5

    # Check for exactly 4 units
    if len(monosaccharide_matches) != 4:
        return False, f"Found {len(monosaccharide_matches)} monosaccharide units, requires 4 for tetrasaccharide"
    
    # Count glycosidic bonds (C-O-C connecting rings)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4]")
    glycosidic_bond_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)

    # A tetrasaccharide will have 3 glycosidic bonds
    if len(glycosidic_bond_matches) < 3 :
         return False, f"Found {len(glycosidic_bond_matches)} glycosidic bonds, requires 3 for tetrasaccharide"

    return True, "Contains 4 monosaccharide units connected via glycosidic bonds"