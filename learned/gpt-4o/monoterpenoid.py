"""
Classifies: CHEBI:25409 monoterpenoid
"""
from rdkit import Chem

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is classified as a monoterpenoid based on its SMILES string.
    Monoterpenoids are derived from monoterpenes, typically having a C10 skeleton that may be rearranged or modified.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Allow a range for C10 to account for rearrangement or slight modifications
    if c_count < 9 or c_count > 12:
        return False, f"Carbon count of {c_count} falls outside typical monoterpenoid range (9-12)"

    # Check for cyclic structures (commonly found in monoterpenoids)
    ring_info = mol.GetRingInfo()
    if not ring_info.IsInitialized():
        return False, "Ring information could not be determined"
    num_rings = len(ring_info.AtomRings())
    if num_rings == 0:
        return False, "No ring structures detected"

    # (Optionally) Check for functional groups typical of monoterpenoids (OH, ketones, etc.)
    alcohol_pattern = Chem.MolFromSmarts("[CX3][OH]")
    ketone_pattern = Chem.MolFromSmarts("C(=O)[C]")
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C]")
    if not (mol.HasSubstructMatch(alcohol_pattern) or mol.HasSubstructMatch(ketone_pattern) or mol.HasSubstructMatch(ester_pattern)):
        return False, "No characteristic functional groups for monoterpenoids found (alcohol, ketone, ester)"

    return True, "Matches C10 skeleton estimate with typical monoterpenoid functional groups"