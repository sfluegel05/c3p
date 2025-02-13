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

    # Check for at least three hydroxyl groups (monosaccharides typically have multiple -OH groups)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 3:
        return False, "Less than 3 hydroxyl groups"

    # Detect cyclic structures and presence of carbonyl
    furanose_pattern = Chem.MolFromSmarts("O1[C@@H]([*])[C@@H]([*])O[C@H]1")
    pyranose_pattern = Chem.MolFromSmarts("O1[C@H]([*])[C@H]([*])O[C@@H]([*])O1")
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")

    if not any(mol.HasSubstructMatch(p) for p in [furanose_pattern, pyranose_pattern, carbonyl_pattern]):
        return False, "No indicative cyclic structure or carbonyl found for monosaccharide"

    # Ensure no glycosidic connections
    anomeric_carbon_pattern = Chem.MolFromSmarts("[C@@H]1O[C@H](OC=*1)[*]")
    if mol.HasSubstructMatch(anomeric_carbon_pattern):
        return False, "Contains glycosidic connections"

    return True, "Matches structure of a monosaccharide with polyhydroxyl groups and cyclic or linear carbonyl group"