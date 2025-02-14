"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: CHEBI:15955 monosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    A monosaccharide is a polyhydroxy aldehyde or polyhydroxy ketone with three or more carbon atoms,
    without glycosidic connection to other units.

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

    # Count carbons, hydrogens, oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Must have at least 3 carbons
    if c_count < 3:
        return False, "Too few carbon atoms for monosaccharide"

    # Check for carbonyl group (aldehyde or ketone)
    carbonyl_pattern = Chem.MolFromSmarts("C(=O)")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl group found, not an aldehyde or ketone"

    # Count hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))

    # Check for polyhydroxy (at least 2 hydroxyls, allowing for deoxy sugars)
    if hydroxyl_count < 2:
        return False, "Not enough hydroxyl groups for monosaccharide"

    # Check for glycosidic bonds (shouldn't have any)
    glycosidic_pattern = Chem.MolFromSmarts("[OX2][CX4][OX2]")
    if mol.HasSubstructMatch(glycosidic_pattern):
        return False, "Contains glycosidic bonds, not a monosaccharide"

    # Check for aldose, ketose, furanose, or pyranose patterns
    aldose_pattern = Chem.MolFromSmarts("[CH2](C(=O))[CH](O)[CH](O)[CH](O)[CH](O)[CH](O)*")
    ketose_pattern = Chem.MolFromSmarts("[CH2]([CH](O)[CH](O))[CH](O)[CH](=O)[CH](O)[CH](O)*")
    furanose_pattern = Chem.MolFromSmarts("O1[CH](O)[CH](O)[CH](O)[CH](O)[CH]1")
    pyranose_pattern = Chem.MolFromSmarts("O1[CH](O)[CH](O)[CH](O)[CH](O)[CH](O)[CH]1")
    if not (mol.HasSubstructMatch(aldose_pattern) or mol.HasSubstructMatch(ketose_pattern) or
            mol.HasSubstructMatch(furanose_pattern) or mol.HasSubstructMatch(pyranose_pattern)):
        return False, "Does not match aldose, ketose, furanose, or pyranose pattern"

    # Additional checks based on counts
    if h_count < 2 * o_count:
        return False, "Too few hydrogens for the number of oxygens"
    if c_count > 8:
        return False, "Carbon chain too long for monosaccharide"

    return True, "Contains carbonyl group and at least 2 hydroxyls without glycosidic bonds, matches monosaccharide pattern"