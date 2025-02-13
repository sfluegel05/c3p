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
    
    # Look for steroid backbone pattern (cyclopenta[a]phenanthrene nucleus)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4C5CCCC(C5)C3CCC12C4")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus found"
    
    # Look for lactone ring pattern (cyclic ester)
    lactone_pattern = Chem.MolFromSmarts("C(=O)O[C@@H]1O")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"
    
    # Look for sugar moieties (cyclic with multiple oxygen atoms)
    sugar_pattern = Chem.MolFromSmarts("C1(CO)O[C@@H]([C@@H]([C@@H]1O)O)CO")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar moieties found"
    
    return True, "Contains steroid nucleus with lactone ring and sugar moieties indicative of cardiac glycoside"