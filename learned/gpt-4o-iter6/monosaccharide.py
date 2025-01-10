"""
Classifies: CHEBI:35381 monosaccharide
"""
from rdkit import Chem

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    Monosaccharides are polyhydroxy aldehydes or ketones with three or more carbon atoms,
    plus they may appear as cyclic hemiacetals or hemiketals.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Monosaccharides must have at least 3 carbon atoms
    if c_count < 3:
        return False, "Too few carbon atoms"

    # Check for hydroxyl groups (monosaccharides usually have several -OH groups)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Less than 2 hydroxyl groups"

    # Detect ring structures (capture furanose, pyranose rings, etc.)
    furanose_pattern = Chem.MolFromSmarts("O1[C@@H]([*])C([*])O[*]1")
    pyranose_pattern = Chem.MolFromSmarts("O1[C@H]([*])[C@H]([*])O[C@@H]([*])O1")
    has_ring = mol.HasSubstructMatch(furanose_pattern) or mol.HasSubstructMatch(pyranose_pattern)

    # Verify presence of either a linear carbonyl group or a ring
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    if not mol.HasSubstructMatch(carbonyl_pattern) and not has_ring:
        return False, "No carbonyl or cyclic structure indicative of monosaccharide"

    # Ensure no glycosidic connections
    # Check for any potential ether linkages involving anomeric carbons
    ether_linkage_pattern = Chem.MolFromSmarts("[OD2]([C@H]1)[*]")
    if mol.HasSubstructMatch(ether_linkage_pattern):
        return False, "Contains glycosidic connection"

    return True, "Matches structure of a monosaccharide with polyhydroxyl groups and potentially hidden cyclic carbonyl group"