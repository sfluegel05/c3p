"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    A diterpenoid has a C20 skeleton derived from a diterpene, potentially rearranged or modified.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # More flexible carbon count check
    if c_count < 16 or c_count > 40:
        return False, f"Uncommon carbon count ({c_count}) for diterpenoids"
    
    # Check for presence of ring structures
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 1:
        return False, "Diterpenoids typically have ring structures"

    # Check for presence of double bonds (C=C, includes alkene SMARTS)
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() == 2)
    if double_bond_count < 1:
        return False, "Few double bonds; uncommon for diterpenoids"

    # Look for functional groups using SMARTS patterns typical for diterpenoids
    patterns = {
        "hydroxyl": Chem.MolFromSmarts("[CX4][OX2H]"),
        "carbonyl": Chem.MolFromSmarts("[CX3](=O)"),
        "ether": Chem.MolFromSmarts("[OX2][CX4]"),
        "epoxide": Chem.MolFromSmarts("C1OC1"),
        "methyl": Chem.MolFromSmarts("[CH3]"),
    }

    for name, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            break
    else:
        return False, "Missing typical functional groups like hydroxyls, epoxides, or methyl groups"

    # Check stereochemistry, including unassigned
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 1:
        return False, "Lack of chiral centers; uncommon for diterpenoids"

    return True, "Molecule matches diterpenoid characteristics with flexible carbon count, ring structures, and typical functional groups"