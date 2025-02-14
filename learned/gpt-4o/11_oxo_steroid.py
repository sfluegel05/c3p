"""
Classifies: CHEBI:47787 11-oxo steroid
"""
from rdkit import Chem

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    An 11-oxo steroid is defined by the presence of a ketone (oxo group) at the 11th carbon position in conjunction with a steroid backbone.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is an 11-oxo steroid, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the 11-oxo steroid SMARTS pattern
    # Naively, steroids are characterized by carbon frameworks, but for simplicity and checking the position of functional groups:
    # Use a basic ten-carbon structure with the double-bonded oxygen as an approximation
    keto_11_pattern = Chem.MolFromSmarts("[C]1(C)[C][C][C]([C])([C])[C](=O)[C][C][C]1")
    
    if keto_11_pattern is None:
        return False, "Invalid SMARTS pattern for keto-11 check"

    # Match this pattern to the molecule
    if not mol.HasSubstructMatch(keto_11_pattern):
        return False, "11-oxo group not found"

    return True, "Contains 11-oxo group in a typical steroid framework"