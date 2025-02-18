"""
Classifies: CHEBI:23044 carotenoid
"""
"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    Carotenoids are tetraterpenoids (C40-derived) with conjugated polyene chains, excluding retinoids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check number of double bonds (at least 8 for conjugated system)
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bonds < 8:
        return False, f"Found {double_bonds} double bonds (needs ≥8)"

    # Check carbon count (derived from C40 backbone, allowing modifications)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30:
        return False, f"Only {c_count} carbons (needs ≥30)"

    # Check for methyl groups adjacent to conjugated systems (isoprenoid feature)
    methyl_pattern = Chem.MolFromSmarts("[CH3]-[C]=,*")  # Methyl near conjugated bonds
    if not mol.HasSubstructMatch(methyl_pattern):
        return False, "No methyl groups near conjugated system"

    # Molecular weight check to exclude retinoids (carotenoids typically >400 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da)"

    # Check for beta-ionone ring (retinoid marker) with short chain exclusion
    beta_ionone = Chem.MolFromSmarts("[C]1([C])=C([CH3])C([C])([C])CCC1")
    if mol.HasSubstructMatch(beta_ionone) and c_count < 35:
        return False, "Beta-ionone ring with short chain (retinoid)"

    return True, "Conjugated polyene chain with isoprenoid features and sufficient size"