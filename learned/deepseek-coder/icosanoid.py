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

    # Check for approximately 20 carbons in the molecule (18-22)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (18 <= c_count <= 22):
        return False, f"Expected 18-22 carbons, found {c_count}"

    # Check for at least one double bond (polyunsaturated)
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bond_count < 1:
        return False, f"Expected at least one double bond, found {double_bond_count}"

    # Check for functional groups: hydroxyl, hydroperoxy, carboxyl, or their derivatives
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroperoxy_pattern = Chem.MolFromSmarts("[OX2][OX1]")
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    ether_pattern = Chem.MolFromSmarts("[OX2]")

    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_hydroperoxy = mol.HasSubstructMatch(hydroperoxy_pattern)
    has_carboxyl = mol.HasSubstructMatch(carboxyl_pattern)
    has_ester = mol.HasSubstructMatch(ester_pattern)
    has_ether = mol.HasSubstructMatch(ether_pattern)

    if not (has_hydroxyl or has_hydroperoxy or has_carboxyl or has_ester or has_ether):
        return False, "Missing characteristic functional groups (hydroxyl, hydroperoxy, carboxyl, ester, or ether)"

    # Check molecular weight (icosanoids typically have MW > 250 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, f"Molecular weight {mol_wt:.2f} Da is too low for an icosanoid"

    # Check for characteristic icosanoid structure: long carbon chain with functional groups
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Missing characteristic long carbon chain structure"

    return True, "Contains characteristic icosanoid structure: 18-22 carbon backbone with functional groups and at least one double bond"