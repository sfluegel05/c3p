"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    Diterpenoids are characterized by a core derived from 20 carbon atoms of a diterpene, often modified or rearranged.

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

    # Count number of carbon atoms; allow flexibility due to potential rearrangements
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15 or c_count > 30:
        return False, "Carbon count not within extended diterpenoid range (15-30)"

    # Retrieve ring information to check for potential diterpenoid structures
    ring_info = mol.GetRingInfo()
    if not ring_info.IsInitialized():
        return False, "No ring data available"
    ring_sizes = [len(ring) for ring in ring_info.AtomRings() if 5 <= len(ring) <= 8]
    if len(ring_sizes) < 1:
        return False, "Lacks typical ring structures of diterpenoids, expected 5-8 members"

    # Check for presence of relevant functional groups (hydroxyl, carbonyl, ether, epoxide)
    functional_group_pattern_1 = Chem.MolFromSmarts("[OX2H,C]=[OX1,O]")
    functional_group_pattern_2 = Chem.MolFromSmarts("[O,N]")
    if not mol.HasSubstructMatch(functional_group_pattern_1) and not mol.HasSubstructMatch(functional_group_pattern_2):
        return False, "Missing expected functional groups (e.g., hydroxyl, carbonyl, ether, epoxide)"

    return True, "Exhibits typical characteristics of diterpenoids (ring structure and functional modifications)"