"""
Classifies: CHEBI:17389 2-monoglyceride
"""
from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride has a glycerol backbone with an acyl substituent at the 2-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern with acyl group at the 2-position
    # Ensure that the first and third OH groups are free, while the middle is esterified
    glycerol_pattern = Chem.MolFromSmarts("O[C@@H](C(=O)[*])CO")
    
    # There must be a match for the pattern that connects the 2-position carbon with an acyl group
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No acyl substituent at the 2-position via ester linkage"
    
    # Check for absence of large chains or additional substituents that don't fit the monoglyceride pattern
    # We will exclude any match with additional phosphates or large complex substituents
    complex_substituents = Chem.MolFromSmarts("P|N|S")
    if mol.HasSubstructMatch(complex_substituents):
        return False, "Contains complex substituents that are not characteristic of a 2-monoglyceride"
    
    return True, "Contains glycerol backbone with an acyl substituent at the 2-position (ester linkage)"