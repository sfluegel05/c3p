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

    # Check for multiple double bonds (typically 2-6)
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bonds < 2 or double_bonds > 6:
        return False, f"Expected 2-6 double bonds, found {double_bonds}"

    # Check for oxygen-containing functional groups (hydroxyl, carboxyl, epoxide)
    oxygen_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_atoms < 2:
        return False, f"Expected at least 2 oxygen atoms, found {oxygen_atoms}"

    # Check molecular weight (should be around 300-400 Da for C20 molecules)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 400:
        return False, f"Molecular weight {mol_wt:.2f} is outside the expected range for icosanoids"

    # Check for specific functional groups (hydroxyl, carboxyl, epoxide)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    epoxide_pattern = Chem.MolFromSmarts("[OX2]1[CX4][CX4]1")

    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_carboxyl = mol.HasSubstructMatch(carboxyl_pattern)
    has_epoxide = mol.HasSubstructMatch(epoxide_pattern)

    if not (has_hydroxyl or has_carboxyl or has_epoxide):
        return False, "No hydroxyl, carboxyl, or epoxide groups found"

    return True, "Contains 20-carbon backbone with multiple double bonds and oxygen-containing functional groups"