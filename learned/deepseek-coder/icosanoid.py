"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: CHEBI:36080 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    An icosanoid is a signaling molecule derived from the oxidation of C20 essential fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 20 carbons in the molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 20:
        return False, f"Expected 20 carbons, found {c_count}"

    # Check for multiple double bonds (polyunsaturated)
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bond_count < 2:
        return False, f"Expected multiple double bonds, found {double_bond_count}"

    # Check for functional groups: hydroxyl, hydroperoxy, or carboxyl
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroperoxy_pattern = Chem.MolFromSmarts("[OX2][OX1]")
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")

    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_hydroperoxy = mol.HasSubstructMatch(hydroperoxy_pattern)
    has_carboxyl = mol.HasSubstructMatch(carboxyl_pattern)

    if not (has_hydroxyl or has_hydroperoxy or has_carboxyl):
        return False, "Missing hydroxyl, hydroperoxy, or carboxyl functional group"

    # Check molecular weight (icosanoids typically have MW around 300-400 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 400:
        return False, f"Molecular weight {mol_wt:.2f} Da is outside the typical range for icosanoids"

    return True, "Contains 20-carbon backbone with multiple double bonds and functional groups (hydroxyl, hydroperoxy, or carboxyl)"