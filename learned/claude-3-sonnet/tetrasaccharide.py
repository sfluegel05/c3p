"""
Classifies: CHEBI:50126 tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetrasaccharide(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide consists of exactly 4 monosaccharide units connected by glycosidic bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple[bool, str]: (is_tetrasaccharide, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # More flexible thresholds
    if c_count < 12:  # Allow for modified tetrasaccharides
        return False, "Too few carbons for a tetrasaccharide"
    if c_count > 40:  # Limit maximum carbons
        return False, "Too many carbons for a tetrasaccharide"
    if o_count < 8:  # Allow for modified tetrasaccharides
        return False, "Too few oxygens for a tetrasaccharide"
    if o_count > 30:  # Limit maximum oxygens
        return False, "Too many oxygens for a tetrasaccharide"

    # Count ring systems
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    
    if ring_count < 3 or ring_count > 6:  # Allow some flexibility but not too much
        return False, f"Found {ring_count} rings, expect around 4 for tetrasaccharide"

    # More specific sugar ring patterns
    pyranose = Chem.MolFromSmarts("[C]1[C][C][C]([O,C])[C]O1")
    furanose = Chem.MolFromSmarts("[C]1[C][C]([O,C])[C]O1")
    
    pyranose_matches = len(mol.GetSubstructMatches(pyranose))
    furanose_matches = len(mol.GetSubstructMatches(furanose))
    total_sugar_rings = pyranose_matches + furanose_matches
    
    if total_sugar_rings < 3:
        return False, f"Found only {total_sugar_rings} sugar rings, need at least 3"
    if total_sugar_rings > 5:
        return False, f"Found {total_sugar_rings} sugar rings, too many for tetrasaccharide"

    # Improved glycosidic bond pattern
    glycosidic = Chem.MolFromSmarts("[C;R][O][C;!$(C=O)]")  # Exclude ester bonds
    glycosidic_matches = len(mol.GetSubstructMatches(glycosidic))
    
    if glycosidic_matches < 2:
        return False, "Too few glycosidic bonds for tetrasaccharide"
    if glycosidic_matches > 5:
        return False, "Too many glycosidic bonds for tetrasaccharide"

    # Count hydroxyl groups (typical for sugars)
    hydroxyl = Chem.MolFromSmarts("[C][OH]")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl))
    
    if hydroxyl_count < 3:
        return False, "Too few hydroxyl groups for tetrasaccharide"
    
    # Allow more aromatic atoms for modified sugars
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms > 6:  # More lenient threshold
        return False, "Too many aromatic atoms for tetrasaccharide"

    # Look for typical sugar carbon patterns
    sugar_carbon = Chem.MolFromSmarts("[C]([O])([O])[C]([O])")
    sugar_matches = len(mol.GetSubstructMatches(sugar_carbon))
    
    if sugar_matches < 2:
        return False, "Missing typical sugar carbon patterns"

    # Check for anomeric carbons
    anomeric = Chem.MolFromSmarts("[C;R]([O])([O])")
    anomeric_matches = len(mol.GetSubstructMatches(anomeric))
    
    if anomeric_matches < 2:
        return False, "Too few anomeric carbons for tetrasaccharide"
    if anomeric_matches > 6:
        return False, "Too many anomeric carbons for tetrasaccharide"

    return True, "Contains approximately 4 monosaccharide units connected by glycosidic bonds"