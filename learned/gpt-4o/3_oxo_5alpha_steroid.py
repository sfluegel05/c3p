"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
from rdkit import Chem

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if a molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more accurate steroid backbone pattern (four rings)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4CCCC(C4)C3CCC2C1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for the presence of a 3-oxo group in the steroid nucleus
    oxo_pattern = Chem.MolFromSmarts("C(=O)[C@H]1CC[C@H]2")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if not oxo_matches:
        return False, "No 3-oxo group found"

    # Ensure 5alpha configuration by identifying the appropriate stereochemistry
    five_alpha_pattern = Chem.MolFromSmarts("C[C@@H]1CCC2")
    if not mol.HasSubstructMatch(five_alpha_pattern):
        return False, "No 5alpha configuration found"

    # Passed all checks
    return True, "Contains a steroid backbone with a 3-oxo group and the 5alpha stereochemical configuration"