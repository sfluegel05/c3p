"""
Classifies: CHEBI:17389 2-monoglyceride
"""
from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride has a glycerol backbone with an acyl substituent at the 2-position (middle carbon).
    
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

    # Add explicit hydrogens to improve substructure matching
    mol = Chem.AddHs(mol)

    # Define the SMARTS pattern for a 2-monoglyceride
    # Looking for the structure HO-CH2-CH(O-C(=O)-R)-CH2OH
    glycerol_pattern = Chem.MolFromSmarts("OCC(O[C:1](=O)[C:2])CO")
    if glycerol_pattern is None:
        return False, "Error in constructing SMARTS pattern"

    # Check for the presence of acyl substituent at 2-position
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No acyl substituent at the 2-position via ester linkage"

    # Confirm there are no additional functional groups like nitrogen, phosphorus, or sulfur
    complex_substituents = Chem.MolFromSmarts("[#7,#15,#16]")
    if mol.HasSubstructMatch(complex_substituents):
        return False, "Contains complex substituents that are not characteristic of a 2-monoglyceride"

    return True, "Contains glycerol backbone with an acyl substituent at the 2-position (ester linkage)"