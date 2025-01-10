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
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:  # Typical tetrasaccharide has at least 20 carbons
        return False, "Too few carbons for a tetrasaccharide"
    if o_count < 15:  # Typical tetrasaccharide has at least 15 oxygens
        return False, "Too few oxygens for a tetrasaccharide"

    # More flexible pattern for sugar rings (both pyranose and furanose)
    pyranose_pattern = Chem.MolFromSmarts("[C]1[C][C][C]([O,C])[C]O1")
    furanose_pattern = Chem.MolFromSmarts("[C]1[C][C]([O,C])[C]O1")
    
    pyranose_matches = len(mol.GetSubstructMatches(pyranose_pattern))
    furanose_matches = len(mol.GetSubstructMatches(furanose_pattern))
    total_sugar_rings = pyranose_matches + furanose_matches
    
    if total_sugar_rings < 4:
        return False, f"Found only {total_sugar_rings} sugar rings, need at least 4"
    
    # Look for glycosidic bonds between sugar rings
    # Pattern matches C-O-C where one C is part of a ring
    glycosidic_pattern = Chem.MolFromSmarts("[C;R][O][C]")
    glycosidic_matches = len(mol.GetSubstructMatches(glycosidic_pattern))
    
    if glycosidic_matches < 3:
        return False, "Too few glycosidic bonds for a tetrasaccharide"

    # Count oxygen atoms bonded to carbon (better than counting explicit hydroxyls)
    c_o_pattern = Chem.MolFromSmarts("[C][O]")
    c_o_count = len(mol.GetSubstructMatches(c_o_pattern))
    
    if c_o_count < 15:  # Tetrasaccharides typically have many C-O bonds
        return False, "Too few C-O bonds for a tetrasaccharide"

    # Check for characteristics that would exclude it from being a tetrasaccharide
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms > 2:  # Allow some aromatic atoms for modified sugars
        return False, "Too many aromatic atoms for a tetrasaccharide"

    # Look for typical sugar carbon patterns
    sugar_carbon_pattern = Chem.MolFromSmarts("[C]([O])([O])[C]([O])")
    if not mol.HasSubstructMatch(sugar_carbon_pattern):
        return False, "Missing typical sugar carbon pattern"

    # Additional check for anomeric carbons (carbons with two oxygens attached)
    anomeric_pattern = Chem.MolFromSmarts("[C;R]([O])([O])")
    anomeric_matches = len(mol.GetSubstructMatches(anomeric_pattern))
    
    if anomeric_matches < 3:  # Should have at least 3 anomeric carbons
        return False, "Too few anomeric carbons for a tetrasaccharide"

    return True, "Contains 4 monosaccharide units connected by glycosidic bonds"