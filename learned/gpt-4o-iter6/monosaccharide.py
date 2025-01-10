"""
Classifies: CHEBI:35381 monosaccharide
"""
from rdkit import Chem

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    Monosaccharides are polyhydroxy aldehydes or ketones with three or more carbon atoms,
    and may appear as cyclic hemiacetals or hemiketals.

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

    # Check for at least three carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return False, "Too few carbon atoms (minimum 3 needed)"
    
    # Look for three or more hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_count < 3:
        return False, f"Less than 3 hydroxyl groups (found {hydroxyl_count})"

    # Look for cyclic hemiacetal/hemiketal structures (furanose and pyranose forms)
    furanose_pattern = Chem.MolFromSmarts("O1[C@H]([C@@H]([C@H]([C@@H]1[O]))[O])[C@H](O)C")
    pyranose_pattern = Chem.MolFromSmarts("O1[C@H]([C@@H]([C@H]([C@H]([C@@H]1[O]))[O])[O])[C@H](O)C")
    
    # Check for linear aldehyde or ketone presence
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    ketone_pattern = Chem.MolFromSmarts("[#6][C](=O)[#6]")

    # Validate cyclic or linear presence of hemiacetal/hemiketal forms
    if not (mol.HasSubstructMatch(furanose_pattern) or mol.HasSubstructMatch(pyranose_pattern) 
            or mol.HasSubstructMatch(aldehyde_pattern) or mol.HasSubstructMatch(ketone_pattern)):
        return False, "No cyclic or linear carbonyl structures"

    # Validate absence of glycosidic connections
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C@H]1O[C@H](OC)[C@H](O)[C@H]([C@H]1O)C")
    if mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "Contains glycosidic connections, indicating non-monosaccharide"

    return True, "Potential monosaccharide based on polyhydroxyl groups and structural elements"