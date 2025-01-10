"""
Classifies: CHEBI:35381 monosaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    Monosaccharides are polyhydroxy aldehydes or ketones, with at least three carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return False, "Too few carbon atoms; monosaccharides must have at least 3"
    
    # Identify the number of hydroxyl groups (O-H patterns)
    oh_pattern = Chem.MolFromSmarts("[OX2H]")  # Hydroxyl group SMARTS pattern
    hydroxyl_matches = mol.GetSubstructMatches(oh_pattern)
    if len(hydroxyl_matches) < 2:
        return False, f"Insufficient hydroxyl groups found; needed at least 2, found {len(hydroxyl_matches)}"

    # Check for potential carbonyl groups or hemiacetal/hemiketal formations
    carbonyl_pattern = Chem.MolFromSmarts("C=O")  # simple carbonyl
    hemiacetal_pattern = Chem.MolFromSmarts("[C;R][O;R]")  # cyclic ether indicating a hemiacetal
    if not (mol.HasSubstructMatch(carbonyl_pattern) or mol.HasSubstructMatch(hemiacetal_pattern)):
        return False, "No potential carbonyl or cyclic oxy-functional group detected."
    
    # Check for appropriate ring sizes indicative of furanose or pyranose (5 or 6-membered rings)
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.BondRings()]
    if not any(size in [5, 6] for size in ring_sizes):
        return False, "No furanose or pyranose ring size detected."

    return True, "Structure is a monosaccharide with potential cyclic or acyclic carbonyl group and sufficient hydroxyls."