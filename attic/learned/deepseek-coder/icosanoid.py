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

    # Check for approximately 20 carbons in the molecule (16-25)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (16 <= c_count <= 25):
        return False, f"Expected 16-25 carbons, found {c_count}"

    # Check for at least one double bond (polyunsaturated)
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bond_count < 1:
        return False, f"Expected at least one double bond, found {double_bond_count}"

    # Check for characteristic functional groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroperoxy_pattern = Chem.MolFromSmarts("[OX2][OX1]")
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    ether_pattern = Chem.MolFromSmarts("[OX2]")
    cyclopentane_pattern = Chem.MolFromSmarts("C1CCCC1")  # Common in prostaglandins

    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_hydroperoxy = mol.HasSubstructMatch(hydroperoxy_pattern)
    has_carboxyl = mol.HasSubstructMatch(carboxyl_pattern)
    has_ester = mol.HasSubstructMatch(ester_pattern)
    has_ether = mol.HasSubstructMatch(ether_pattern)
    has_cyclopentane = mol.HasSubstructMatch(cyclopentane_pattern)

    if not (has_hydroxyl or has_hydroperoxy or has_carboxyl or has_ester or has_ether or has_cyclopentane):
        return False, "Missing characteristic functional groups or structures"

    # Check molecular weight (icosanoids typically have MW > 250 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, f"Molecular weight {mol_wt:.2f} Da is too low for an icosanoid"

    # Check for characteristic icosanoid structure patterns
    # Pattern 1: Long carbon chain with functional groups
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    # Pattern 2: Cyclopentane ring with adjacent functional groups
    prostaglandin_pattern = Chem.MolFromSmarts("[CX4]1[CX4][CX4][CX4][CX4]1")
    
    if not (mol.HasSubstructMatch(chain_pattern) or mol.HasSubstructMatch(prostaglandin_pattern)):
        return False, "Missing characteristic icosanoid structural patterns"

    return True, "Contains characteristic icosanoid structure: 16-25 carbon backbone with functional groups and at least one double bond"