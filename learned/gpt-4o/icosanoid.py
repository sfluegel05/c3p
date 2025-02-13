"""
Classifies: CHEBI:23899 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    Icosanoids are signaling molecules derived from oxidation of C20 essential fatty acids.

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

    # Check the number of carbon atoms (focus on the C20 backbone)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (18 <= c_count <= 24):  # Extend range slightly for connected extensions
        return False, f"Carbon count ({c_count}) not within expected range for C20 backbone"

    # Check for presence of oxygen atoms (indicative of oxidation patterns)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, f"Insufficient number of oxygens ({o_count}), expect oxidized features in icosanoids"

    # Check for presence of double bonds (characteristic of EFAs)
    double_bonds = len([bond for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() == 2])
    if double_bonds < 4:
        return False, f"Insufficient number of double bonds ({double_bonds}) for an icosanoid"

    # Check for typical functional groups: keto (C=O) and hydroxyl (OH)
    keto_pattern = Chem.MolFromSmarts("C=O")
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    if not (mol.HasSubstructMatch(keto_pattern) or mol.HasSubstructMatch(hydroxyl_pattern)):
        return False, "Missing common functional groups (keto/hydroxyl) in icosanoids"

    # Check for possible cyclic structures like cyclopentane rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 1:
        return False, "Cyclopentane or similar ring structure not found in icosanoids"

    # Consider additional linkages (e.g., esterified forms)
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    if mol.HasSubstructMatch(ester_pattern) or mol.HasSubstructMatch(amide_pattern):
        return True, "Structure contains icosanoid features with ester or amide linkages"

    return True, "Molecule generally fits structural criteria suggestive of an icosanoid"