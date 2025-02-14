"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    A diterpenoid typical core has 20 carbons derived from a C20 skeleton,
    modified with rearranged or functional groups.

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
    # Allow more flexibility in carbon count due to modifications
    if c_count < 15 or c_count > 25:
        return False, "Carbon count not within typical diterpenoid range (15-25)"

    # Check for presence of carbon rings
    ring_info = mol.GetRingInfo()
    if not ring_info.IsInitialized():
        return False, "No ring data available"
    ring_sizes = [len(ring) for ring in ring_info.AtomRings() if 5 <= len(ring) <= 7]
    if len(ring_sizes) < 1:
        return False, "Lacks typical ring structures of diterpenoids"

    # Broaden detection for functional groups or modifications (hydroxyl, oxo, epoxide)
    functional_group_pattern = Chem.MolFromSmarts("[!#6][OH0,O]")
    if not mol.HasSubstructMatch(functional_group_pattern):
        return False, "Missing functional groups typical in diterpenoids"

    return True, "Exhibits typical characteristics of diterpenoids (rings and functional groups)"