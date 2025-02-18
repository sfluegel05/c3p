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

    # Check for core icosanoid structure (carbon count around 20)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18 or c_count > 22:
        return False, f"Expected around 20 carbons, found {c_count}"

    # Check for multiple double bonds (at least 2)
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bonds < 2:
        return False, f"Expected at least 2 double bonds, found {double_bonds}"

    # Check for oxygen-containing functional groups (at least 2)
    oxygen_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_atoms < 2:
        return False, f"Expected at least 2 oxygen atoms, found {oxygen_atoms}"

    # More flexible patterns for characteristic icosanoid structures
    # 1. Cyclopentane ring with oxygen-containing groups (prostaglandin-like)
    prostaglandin_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C]1([C]=O)([OH])")
    # 2. Conjugated polyene system (leukotriene-like)
    leukotriene_pattern = Chem.MolFromSmarts("[C]=[C][C]=[C][C]=[C]")
    # 3. Epoxide structure (common in EETs)
    epoxide_pattern = Chem.MolFromSmarts("[OX2]1[CX4][CX4]1")
    # 4. Hydroperoxide structure (common in HPETEs)
    hydroperoxide_pattern = Chem.MolFromSmarts("[OH][OX2]")
    # 5. Carboxyl group (common in many icosanoids)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")

    has_prostaglandin = mol.HasSubstructMatch(prostaglandin_pattern)
    has_leukotriene = mol.HasSubstructMatch(leukotriene_pattern)
    has_epoxide = mol.HasSubstructMatch(epoxide_pattern)
    has_hydroperoxide = mol.HasSubstructMatch(hydroperoxide_pattern)
    has_carboxyl = mol.HasSubstructMatch(carboxyl_pattern)

    if not (has_prostaglandin or has_leukotriene or has_epoxide or has_hydroperoxide or has_carboxyl):
        return False, "No characteristic icosanoid structural patterns found"

    # Check for specific functional groups (hydroxyl, carboxyl, epoxide, peroxide, ester, ether)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")
    ether_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4]")

    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_ester = mol.HasSubstructMatch(ester_pattern)
    has_ether = mol.HasSubstructMatch(ether_pattern)

    if not (has_hydroxyl or has_carboxyl or has_ester or has_ether):
        return False, "No hydroxyl, carboxyl, ester, or ether groups found"

    return True, "Contains characteristic icosanoid structure with multiple double bonds and oxygen-containing functional groups"