"""
Classifies: CHEBI:83970 cardiac glycoside
"""
from rdkit import Chem

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    A cardiac glycoside is characterized by a steroid core, lactone ring, and sugar moieties.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiac glycoside, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for an extended steroid backbone pattern that includes achievable variants
    steroid_pattern = Chem.MolFromSmarts("C1CCC2(C(C1)CCC3C2CCC4(C3CCC4)C5CC=CC=C5)")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus found"

    # Look for common lactone ring patterns
    lactone_pattern_5 = Chem.MolFromSmarts("O=C1OC2C(=O)CCC12")  # 5-membered cyclic ester
    lactone_pattern_6 = Chem.MolFromSmarts("O=C1OC2C(=O)CCCC12")  # 6-membered cyclic ester
    if not (mol.HasSubstructMatch(lactone_pattern_5) or mol.HasSubstructMatch(lactone_pattern_6)):
        return False, "No lactone ring found"

    # Look for broad sugar patterns to capture different sugar ring sizes
    sugar_pattern_pyran = Chem.MolFromSmarts("C1OC([C@@H]([C@H](O1)C)O)C")  # Pyranose conditional pattern
    sugar_pattern_furan = Chem.MolFromSmarts("C1OC([C@@H](O)[C@@H](O1)C)O")  # Furanose conditional pattern
    if not (mol.HasSubstructMatch(sugar_pattern_pyran) or mol.HasSubstructMatch(sugar_pattern_furan)):
        return False, "No sugar moieties found"

    return True, "Contains steroid nucleus with lactone ring and sugar moieties indicative of cardiac glycoside"