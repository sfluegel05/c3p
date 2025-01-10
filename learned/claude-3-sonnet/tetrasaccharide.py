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

    # Basic checks for carbohydrate characteristics
    # Count C, O, H atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:  # Typical tetrasaccharide has at least 24 carbons
        return False, "Too few carbons for a tetrasaccharide"
    if o_count < 15:  # Typical tetrasaccharide has at least 20 oxygens
        return False, "Too few oxygens for a tetrasaccharide"

    # Look for pyranose/furanose rings
    sugar_ring_pattern = Chem.MolFromSmarts("[C]1[C][C][C]([O,C])[C]([O,C])O1")
    ring_matches = mol.GetSubstructMatches(sugar_ring_pattern)
    
    if len(ring_matches) < 4:
        return False, f"Found only {len(ring_matches)} sugar rings, need exactly 4"
    
    # Look for glycosidic bonds (C-O-C between rings)
    glycosidic_pattern = Chem.MolFromSmarts("[C][O][C]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    
    if len(glycosidic_matches) < 3:
        return False, "Too few glycosidic bonds for a tetrasaccharide"

    # Look for hydroxyl groups (typical for sugars)
    hydroxyl_pattern = Chem.MolFromSmarts("[O][H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    if len(hydroxyl_matches) < 8:  # Tetrasaccharides typically have many OH groups
        return False, "Too few hydroxyl groups for a tetrasaccharide"

    # Count ring systems using RDKit's built-in function
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    
    if ring_count < 4:
        return False, f"Found only {ring_count} rings, need exactly 4"
    elif ring_count > 6:  # Allow some flexibility for fused rings
        return False, f"Too many rings ({ring_count}) for a tetrasaccharide"

    # Check for characteristics that would exclude it from being a tetrasaccharide
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms > 2:  # Allow some aromatic atoms for modified sugars
        return False, "Too many aromatic atoms for a tetrasaccharide"

    # Additional check for typical sugar hydroxyl pattern
    sugar_hydroxyl_pattern = Chem.MolFromSmarts("[C]([O][H])([O][H])[C]([O][H])")
    if not mol.HasSubstructMatch(sugar_hydroxyl_pattern):
        return False, "Missing typical sugar hydroxyl pattern"

    # If all checks pass, it's likely a tetrasaccharide
    return True, "Contains 4 monosaccharide units connected by glycosidic bonds"