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
    if c_count < 18 or c_count > 30:  # Allow flexibility for known diterpenoid variants and rearrangements
        return False, f"Uncommon carbon count ({c_count}) for diterpenoids"

    # Check for presence of multiple cycles (ring structures)
    ring_info = mol.GetRingInfo()
    if not ring_info or ring_info.NumRings() < 2:
        return False, "Diterpenoids typically have multiple ring structures"
    
    # Check for presence of double bonds (C=C), common in terpenoid structures
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() == 2)
    if double_bond_count < 2:
        return False, "Few double bonds; uncommon for diterpenoids"

    # Check for functional groups like alcohols, ketones, or ethers typical in diterpenoids
    has_oxygen = any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())
    if not has_oxygen:
        return False, "Lack of typical functional groups like alcohols or ethers"

    # Check for common diterpenoid patterns, such as epoxides or specific hydroxyl arrangements
    diterpenoid_patterns = [
        Chem.MolFromSmarts("[C@](C)(O)"),  # Example: Hydroxyl group pattern
        Chem.MolFromSmarts("O=[C]OC"),   # Example: Ester or epoxide pattern
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in diterpenoid_patterns):
        return False, "Missing signature substructures for diterpenoids"

    return True, "Molecule matches diterpenoid characteristics with flexible carbon count and typical structural features"