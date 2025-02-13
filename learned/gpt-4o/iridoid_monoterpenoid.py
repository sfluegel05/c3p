"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    Iridoid monoterpenoids often include a cyclopentane ring fused to a six-membered oxygen heterocycle.

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

    # Pattern for cyclopentane ring
    cyclopentane_pattern = Chem.MolFromSmarts("[CH2]1[CH2][CH2][CH2][CH2]1")
    if not mol.HasSubstructMatch(cyclopentane_pattern):
        return False, "No cyclopentane ring found"
    
    # Pattern for six-membered oxygen heterocycle (pyran or furan could be alternatives)
    # Use a wildcard pattern with an oxygen in it
    oxygen_heterocycle_pattern = Chem.MolFromSmarts("[O;R]1[CH2][CH2][CH2][CH2][CH2]1")
    if not mol.HasSubstructMatch(oxygen_heterocycle_pattern):
        return False, "No six-membered oxygen heterocycle found"

    # Pattern for ring fusion: Physical overlap in structure typically with shared atoms
    # For simplicity: combining patterns
    fused_ring_pattern = Chem.MolFromSmarts("C1C(C2)CCC2O1")
    if not mol.HasSubstructMatch(fused_ring_pattern):
        return False, "No fusion between cyclopentane and six-membered oxygen heterocycle found"

    return True, "Structure consists of a cyclopentane ring fused to a six-membered oxygen heterocycle, characteristic of iridoid monoterpenoids"