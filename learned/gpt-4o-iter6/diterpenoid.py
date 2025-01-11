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

    # Broaden carbon count range based on diterpenoids' variability
    if c_count < 15 or c_count > 50:
        return False, f"Uncommon carbon count ({c_count}) for diterpenoids"
    
    # Allow for diverse structures including cyclic, bicyclic, etc.
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 1 and c_count < 20:
        return False, "Structures like cyclic or polycyclic are typical for diterpenoids"

    # Check for double bonds
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() == 2)
    if double_bond_count < 1:  # May check minimum double bond occurrence for diterpenoids
        return False, "Few double bonds; uncommon for diterpenoids"

    # Look for common functional groups
    possible_groups = [
        "[CX4][OX2H]",  # hydroxyl
        "[CX3](=O)",    # carbonyl
        "[OX2][CX4]",   # ether
        "C1OC1",        # epoxide
        "[CH3]",        # methyl
        "[OX2H]",       # hydroxyl group
    ]

    if not any(mol.HasSubstructMatch(Chem.MolFromSmarts(p)) for p in possible_groups):
        return False, "Missing typical functional groups like hydroxyls, carbonyls, ethers, epoxides, or methyl groups"

    # Ensure presence of at least one chiral center
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 1:
        return False, "Lack of chiral centers; uncommon for diterpenoids"

    return True, "Molecule matches diterpenoid characteristics with flexible carbon count, diverse structures, and typical functional groups"