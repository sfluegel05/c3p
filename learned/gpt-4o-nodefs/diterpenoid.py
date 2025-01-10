"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    Diterpenoids typically contain 20 carbons (allowing for some variation) and exhibit complex ring structures.

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

    # Count carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 18 or num_carbons > 22:
        return False, f"Contains {num_carbons} carbon atoms, typical diterpenoids vary around 20"

    # Get the ring information about the molecule
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 2:
        return False, "Diterpenoids typically have multiple ring structures"

    # Check for presence of diverse functional groups
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
    ketone_pattern = Chem.MolFromSmarts("C=O")
    ether_pattern = Chem.MolFromSmarts("[OX2]([#6])[#6]")
    
    has_alcohol = mol.HasSubstructMatch(alcohol_pattern)
    has_ketone = mol.HasSubstructMatch(ketone_pattern)
    has_ether = mol.HasSubstructMatch(ether_pattern)

    if not (has_alcohol or has_ketone or has_ether):
        return False, "Lacks common functional groups like alcohols, ketones, or ethers"

    return True, "Contains characteristics typical of diterpenoids, such as complex ring structures and varied functional groups"