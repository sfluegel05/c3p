"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    A diterpenoid typically has a core terpenoid structure with 20 carbons,
    potentially rearranged or functionalized.

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

    # Count number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18 or c_count > 22:
        return False, "Carbon count not within typical diterpenoid range (18-22)"

    # Look for cyclic structures typical in terpenoids
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if len(ring_sizes) < 2 or not any(size >= 5 for size in ring_sizes):
        return False, "Lacks ring structures typical of diterpenoids"

    # Check for specific functional groups like hydroxyl groups and epoxides
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl groups found, typical in diverse diterpenoids"

    # The presence of a five or six-membered carbon ring along with complex structure
    complex_ring_pattern = Chem.MolFromSmarts("C1CCCC1 | C1CCCCC1")
    if not mol.HasSubstructMatch(complex_ring_pattern):
        return False, "Missing complex carbon ring patterns typical for this class"

    return True, "Exhibits characteristics of diterpenoids (C20 skeleton, rings, and functional groups)"