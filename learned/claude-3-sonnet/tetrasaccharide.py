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

    # Basic size check
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:  # Minimum for 4 sugar units
        return False, "Too few carbons for a tetrasaccharide"
    if o_count < 10:  # Minimum oxygens for 4 sugar units
        return False, "Too few oxygens for a tetrasaccharide"

    # Define sugar ring patterns more precisely
    pyranose = Chem.MolFromSmarts("[C]1[C][C][C]([O,N])[C]O1")  # Allow N for amino sugars
    furanose = Chem.MolFromSmarts("[C]1[C][C]([O,N])[C]O1")
    
    pyranose_matches = len(mol.GetSubstructMatches(pyranose))
    furanose_matches = len(mol.GetSubstructMatches(furanose))
    total_sugar_rings = pyranose_matches + furanose_matches
    
    if total_sugar_rings != 4:
        return False, f"Found {total_sugar_rings} sugar rings, need exactly 4 for tetrasaccharide"

    # Check glycosidic linkages
    glycosidic = Chem.MolFromSmarts("[C;R]1[O,N][C;R]")  # Ring carbon to oxygen to ring carbon
    glycosidic_matches = len(mol.GetSubstructMatches(glycosidic))
    
    if glycosidic_matches != 3:  # Exactly 3 linkages needed to connect 4 units
        return False, f"Found {glycosidic_matches} glycosidic bonds, need exactly 3"

    # Verify sugar characteristics
    # Look for characteristic hydroxyl pattern on sugar rings
    sugar_hydroxyl = Chem.MolFromSmarts("[C;R]([O;!R])([O,N;R])")
    hydroxyl_matches = len(mol.GetSubstructMatches(sugar_hydroxyl))
    
    if hydroxyl_matches < 8:  # Expect at least 2 hydroxyls per sugar unit
        return False, "Insufficient hydroxyl groups for tetrasaccharide"

    # Check for anomeric carbons
    anomeric = Chem.MolFromSmarts("[C;R]1[O,N][C;!R]-[O,N]")
    anomeric_matches = len(mol.GetSubstructMatches(anomeric))
    
    if anomeric_matches < 3:  # Need at least 3 anomeric centers
        return False, "Too few anomeric carbons"

    # Verify carbon chain pattern typical of sugars
    sugar_carbon_chain = Chem.MolFromSmarts("[C;R]-[C;R]-[C;R]([O,N])-[C;R]")
    chain_matches = len(mol.GetSubstructMatches(sugar_carbon_chain))
    
    if chain_matches < 4:  # Need at least one per sugar unit
        return False, "Missing characteristic sugar carbon patterns"

    # Check for non-sugar aromatic systems
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms > 2:  # Allow minimal aromatic modification
        return False, "Contains too many aromatic atoms for tetrasaccharide"

    return True, "Contains exactly 4 monosaccharide units connected by glycosidic bonds"