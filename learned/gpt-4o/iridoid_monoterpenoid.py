"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    An iridoid monoterpenoid typically contains a cyclopentane ring fused to a six-membered oxygen heterocycle.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an iridoid monoterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for cyclopentane pattern
    cyclopentane_pattern = Chem.MolFromSmarts("C1CCCC1")
    if not mol.HasSubstructMatch(cyclopentane_pattern):
        return False, "No cyclopentane ring found"
    
    # Look for six-membered oxygen heterocycle
    oxygen_heterocycle_pattern = Chem.MolFromSmarts("O1CCCCC1")
    if not mol.HasSubstructMatch(oxygen_heterocycle_pattern):
        return False, "No six-membered oxygen heterocycle found"

    # Check for fusion between cyclopentane and six-membered oxygen heterocycle
    fused_ring_pattern = Chem.MolFromSmarts("C1(C2CCC1)OCCCC2")
    if not mol.HasSubstructMatch(fused_ring_pattern):
        return False, "No fusion between cyclopentane and oxygen heterocycle found"

    # Optional: look for additional group patterns unique to iridoids 
    # (e.g., lactones or specific glycol groups if defined)
    # For now, we keep it general based on the provided definition and examples

    return True, "Structure consists of a cyclopentane ring fused to a six-membered oxygen heterocycle, characteristic of iridoid monoterpenoids"