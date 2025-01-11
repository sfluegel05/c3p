"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: CHEBI:35381 monosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    A monosaccharide is a polyhydroxy aldehyde or ketone with at least three carbon atoms,
    and it should not be part of a larger oligosaccharide or polysaccharide structure.

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
        return False, "Less than three carbon atoms, not a monosaccharide"

    # Check for potential carbonyl group (aldehyde or ketone)
    # Including both explicit and potential (hydrated) forms
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    potential_carbonyl = Chem.MolFromSmarts("[CX4]([OH])[OH]")  # Hydrated carbonyl
    if not (mol.HasSubstructMatch(carbonyl_pattern) or mol.HasSubstructMatch(potential_carbonyl)):
        return False, "No carbonyl group (aldehyde or ketone) found"

    # Check for hydroxyl groups (at least one required)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1)
    if hydroxyl_count < 1:
        return False, "Insufficient hydroxyl groups for a monosaccharide"

    # Check for glycosidic bonds (more reliable pattern)
    # Looking for O-C-O-C pattern where both oxygens are connected to carbons
    glycosidic_pattern = Chem.MolFromSmarts("[OX2;H0][CX4][CX4][OX2;H0]")
    if mol.HasSubstructMatch(glycosidic_pattern):
        return False, "Glycosidic bond detected, likely part of an oligo/polysaccharide"

    # Check molecular weight (monosaccharides typically < 1000 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 1000:
        return False, "Molecular weight too high for a monosaccharide"

    # Check for typical monosaccharide composition
    # Should have more oxygens than nitrogens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count > o_count:
        return False, "Too many nitrogen atoms for a typical monosaccharide"

    return True, "Polyhydroxy aldehyde or ketone with at least three carbon atoms, no glycosidic bonds, and appropriate composition"