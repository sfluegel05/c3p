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

    # More refined steroid backbone pattern, four rings (A, B, C, D)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4CCCCC4C3CCC2C1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No refined steroid backbone found"

    # Check for 3-oxo group
    oxo_pattern = Chem.MolFromSmarts("C1(=O)[C@@H]2CCC3")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if not oxo_matches:
        return False, "No 3-oxo group found"

    # Ensure 5alpha configuration: A hydrogen atom at the alpha position on C5
    five_alpha_pattern = Chem.MolFromSmarts("C[C@@H]1CCC2")
    if not mol.HasSubstructMatch(five_alpha_pattern):
        return False, "No 5alpha configuration found"

    # Passed all checks
    return True, "Contains a 3-oxo group and the 5alpha stereochemical configuration, indicating a 3-oxo-5alpha-steroid"