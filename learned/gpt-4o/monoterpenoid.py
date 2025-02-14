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

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8 or c_count > 12:
        return False, f"Carbon count of {c_count} falls outside typical monoterpenoid range (8-12)"

    # Check for common monoterpenoid functional groups
    alcohol_pattern = Chem.MolFromSmarts("[CX3][OH]")
    ketone_pattern = Chem.MolFromSmarts("C(=O)[C]")
    ether_pattern = Chem.MolFromSmarts("C-O-C")
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C]")
    if not (mol.HasSubstructMatch(alcohol_pattern) or
            mol.HasSubstructMatch(ketone_pattern) or
            mol.HasSubstructMatch(ether_pattern) or
            mol.HasSubstructMatch(ester_pattern)):
        return False, "No characteristic functional groups for monoterpenoids found (alcohol, ketone, ether, ester)"

    # Check for presence of at least one ring structure
    ring_info = mol.GetRingInfo()
    num_rings = len(ring_info.AtomRings())
    if num_rings > 0:
        return True, "Contains a C10 skeleton with typical monoterpenoid functional groups and ring structure."

    return True, "Contains a C10 skeleton with typical monoterpenoid functional groups."