"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

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
    
    # Diterpenoids typically have a backbone ranging around 20 carbons, but may vary due to rearrangements
    if c_count < 15 or c_count > 25:  # Allow more flexibility than a strict C20
        return False, f"Uncommon carbon count ({c_count}) for diterpenoids"

    # Terpenoid characteristics: made from isoprene (C5) units, but can be rearranged
    # Terpenoids often have C=C bonds, epoxide groups, and can form complex fused rings

    # Check for presence of at least one cycle (ring structure)
    ring_info = mol.GetRingInfo()
    if not ring_info or ring_info.NumRings() < 1:
        return False, "Diterpenoids typically have at least one ring structure"
    
    # Check for presence of double bonds (C=C), common in terpenoid structures
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() == 2)
    if double_bond_count < 2:
        return False, "Few double bonds; uncommon for diterpenoids"

    # Check for functional groups like alcohols, ketones, or ethers typical in diterpenoids
    has_oxygen = any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())
    if not has_oxygen:
        return False, "Lack of typical functional groups like alcohols or ethers"

    return True, "Molecule matches diterpenoid characteristics with flexible carbon count and typical structural features"