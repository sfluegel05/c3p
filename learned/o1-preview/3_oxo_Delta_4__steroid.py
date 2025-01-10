"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
"""
Classifies: 3-oxo-Delta(4) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    A 3-oxo-Delta(4) steroid is defined as a steroid with a ketone group at position 3
    and a double bond between carbons 4 and 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(4) steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid backbone pattern (cyclopentanoperhydrophenanthrene nucleus)
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CCC4C3CCCC4')
    if steroid_pattern is None:
        return False, "Error in steroid backbone SMARTS pattern"
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Define pattern for ketone at position 3 in ring A
    ketone_pattern = Chem.MolFromSmarts('C1=CC(=O)CC[C@H]1*')
    if ketone_pattern is None:
        return False, "Error in ketone SMARTS pattern"
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group at position 3 in ring A found"

    # Define pattern for double bond between carbons 4 and 5 in ring A
    double_bond_pattern = Chem.MolFromSmarts('C1=CC(=O)C=CC1')
    if double_bond_pattern is None:
        return False, "Error in double bond SMARTS pattern"
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No double bond between carbons 4 and 5 in ring A found"

    return True, "Contains steroid backbone with ketone at position 3 and Δ⁴ double bond between carbons 4 and 5"