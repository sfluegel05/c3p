"""
Classifies: CHEBI:35381 monosaccharide
"""
from rdkit import Chem

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    Monosaccharides are polyhydroxy aldehydes or ketones with three or more carbon
    atoms, typically with multiple hydroxyl groups and potentially in cyclic form.

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

    # Count carbon atoms (C)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return False, "Too few carbon atoms for a monosaccharide"

    # Count hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Too few hydroxyl groups for a monosaccharide"

    # Consider the structure as potentially cyclic, check for ring
    ring_count = mol.GetRingInfo().NumRings()
    
    # Allow structures with either explicit carbonyl groups or cyclic forms
    # i.e., hemiacetal/hemiketal should be considered
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=O")
    if mol.HasSubstructMatch(carbonyl_pattern) or ring_count > 0:
        return True, "Structure has characteristic features of monosaccharides (polyhydroxy and cyclic or carbonyl)"

    return False, "Does not match monosaccharide structural requirements"