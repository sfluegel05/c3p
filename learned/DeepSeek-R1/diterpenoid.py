"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: CHEBI:23899 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

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

    # Check carbon count (allowing for modifications)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (18 <= c_count <= 22):
        return False, f"Carbon count {c_count} outside diterpenoid range (18-22)"

    # Check for structural features common in terpenoids
    # 1. Number of rings and double bonds
    rings = mol.GetRingInfo().NumRings()
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if rings + double_bonds < 4:
        return False, f"Insufficient unsaturation (rings: {rings}, DB: {double_bonds})"

    # 2. Presence of methyl groups (common in terpene branching)
    methyl_matches = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CH3]")))
    if methyl_matches < 2:
        return False, f"Only {methyl_matches} methyl groups, expected ≥2"

    # 3. Check for oxygen atoms (common in modified terpenoids)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count == 0:
        return False, "No oxygen atoms present (unlikely modified diterpenoid)"

    return True, f"C{c_count} with {rings} rings, {double_bonds} double bonds, {methyl_matches} methyl groups"