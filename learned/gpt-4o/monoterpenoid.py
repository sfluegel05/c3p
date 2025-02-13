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
    if c_count < 8 or c_count > 20:
        return False, f"Carbon count of {c_count} falls outside typical range for rearranged monoterpenoids (8-20)"

    # Check for functional groups typical of monoterpenoids
    functional_groups = [
        Chem.MolFromSmarts("[CX3][OH]"),     # Alcohol
        Chem.MolFromSmarts("C(=O)[CX3]"),   # Ketone
        Chem.MolFromSmarts("[#6]-O-[#6]"),  # Ether
        Chem.MolFromSmarts("C(=O)O"),       # Acid or Ester
        Chem.MolFromSmarts("C=O"),          # Aldehyde
        Chem.MolFromSmarts("C=C")           # Alkenes indicative of isoprene units
    ]

    if not any(mol.HasSubstructMatch(pattern) for pattern in functional_groups):
        return False, "No characteristic functional groups for monoterpenoids found (alcohol, ketone, ether, ester, acid, aldehyde)"

    # Identify patterns for monoterpenoid skeletons
    terpenoid_patterns = [
        Chem.MolFromSmarts("C1CCC(C)CC1"), # Cyclohexane with branching
        Chem.MolFromSmarts("C=C(C)CC"),    # Isoprenoid unit
        Chem.MolFromSmarts("C1=C(C)CCCC1") # Monoterpene-like carbon rings
    ]

    if not any(mol.HasSubstructMatch(pattern) for pattern in terpenoid_patterns):
        return False, "No typical monoterpenoid carbon skeleton found (rearranged or derived structures)"

    return True, "Contains a rearranged C10 skeleton or functional group characteristic of monoterpenoids."