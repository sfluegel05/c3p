"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: CHEBI:23899 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    Diterpenoids are derived from a C20 diterpene skeleton, which may be modified.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Adjusted carbon count range (allowing for some modification)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (18 <= c_count <= 24):
        return False, f"Carbon count {c_count} outside diterpenoid range (18-24)"

    # Check for oxygen in specific functional groups (hydroxyl, carbonyl, ester)
    oxygen_groups = [
        Chem.MolFromSmarts("[OH]"),  # Hydroxyl
        Chem.MolFromSmarts("C=O"),   # Carbonyl
        Chem.MolFromSmarts("[O][C](=O)"), # Ester or acid
    ]
    has_oxygen = any(mol.HasSubstructMatch(patt) for patt in oxygen_groups if patt)
    if not has_oxygen:
        return False, "No oxygen in key functional groups (hydroxyl/carbonyl/ester)"

    # Check for terpene-like branching patterns (modified from isoprene units)
    # Pattern for branched carbons (at least two branches)
    branch_pattern = Chem.MolFromSmarts("[CH2][C]([CH3])([CH3])[CH2]")
    complex_branch = mol.HasSubstructMatch(branch_pattern)
    
    # Alternative pattern for cyclic terpene structures
    ring_pattern = Chem.MolFromSmarts("[C]1~[C]~[C]~[C]~[C]~[C]1")  # Generic 6-membered ring
    if not (complex_branch or mol.GetRingInfo().NumRings() >= 2):
        return False, "Lacks terpene-like branching or ring structures"

    # Check methyl groups (minimum 1 allowed for modified structures)
    methyl_matches = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CH3]")))
    if methyl_matches < 1:
        return False, f"Only {methyl_matches} methyl groups, expected ≥1"

    # Check unsaturation (combination of rings and double bonds)
    rings = mol.GetRingInfo().NumRings()
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if (rings + double_bonds) < 3:
        return False, f"Insufficient unsaturation (rings: {rings}, DB: {double_bonds})"

    return True, f"C{c_count} with {rings} rings, {double_bonds} double bonds, {methyl_matches} methyl groups"