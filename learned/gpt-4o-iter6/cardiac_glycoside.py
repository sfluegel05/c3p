"""
Classifies: CHEBI:83970 cardiac glycoside
"""
from rdkit import Chem

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    Cardiac glycosides are characterized by a steroid core, lactone ring, and sugar moieties.

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

    # Look for a broad steroid backbone pattern
    # This pattern aims to capture variations in steroid cores found in cardiac glycosides
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3(CCCC4)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus found"

    # Look for common lactone ring patterns (focusing on 5-membered butenolide)
    lactone_pattern = Chem.MolFromSmarts("C1=CC=CC(=O)C1")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"

    # Look for sugar moieties (flexible detection for pyranose and furanose)
    sugar_pattern = Chem.MolFromSmarts("C1OC(C(O)C(O)C1O)C")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar moieties found"

    return True, "Contains steroid nucleus with lactone ring and sugar moieties indicative of cardiac glycoside"