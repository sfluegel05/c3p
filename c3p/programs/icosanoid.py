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

    # Check for 18-22 carbons in the molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18 or c_count > 22:
        return False, f"Expected 18-22 carbons, found {c_count}"

    # Check for multiple double bonds (typically 2-8)
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bonds < 2 or double_bonds > 8:
        return False, f"Expected 2-8 double bonds, found {double_bonds}"

    # Check for oxygen-containing functional groups (hydroxyl, carboxyl, epoxide, peroxide, ester)
    oxygen_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_atoms < 2:
        return False, f"Expected at least 2 oxygen atoms, found {oxygen_atoms}"

    # Check molecular weight (should be around 250-500 Da for C18-C22 molecules)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 500:
        return False, f"Molecular weight {mol_wt:.2f} is outside the expected range for icosanoids"

    # Check for specific functional groups (hydroxyl, carboxyl, epoxide, peroxide, ester)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    epoxide_pattern = Chem.MolFromSmarts("[OX2]1[CX4][CX4]1")
    peroxide_pattern = Chem.MolFromSmarts("[OX2][OX2]")
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")

    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_carboxyl = mol.HasSubstructMatch(carboxyl_pattern)
    has_epoxide = mol.HasSubstructMatch(epoxide_pattern)
    has_peroxide = mol.HasSubstructMatch(peroxide_pattern)
    has_ester = mol.HasSubstructMatch(ester_pattern)

    if not (has_hydroxyl or has_carboxyl or has_epoxide or has_peroxide or has_ester):
        return False, "No hydroxyl, carboxyl, epoxide, peroxide, or ester groups found"

    return True, "Contains 18-22 carbon backbone with multiple double bonds and oxygen-containing functional groups"