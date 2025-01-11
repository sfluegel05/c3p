"""
Classifies: CHEBI:71543 rotenoid
"""
from rdkit import Chem

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    Rotenoids are characterized by a polycyclic structure often with chromene or isoflavone-like backbones with methoxy groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a rotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for rotenoid substructures
    # A proposed pattern might include a chromene unit and a methoxy substituent
    chromene_pattern = Chem.MolFromSmarts("c1ccc2occc2c1")
    methoxy_pattern = Chem.MolFromSmarts("CO")

    # Check for the presence of a chromene-like backbone
    if not mol.HasSubstructMatch(chromene_pattern):
        return False, "No chromene-like backbone found"

    # Check for methoxy substituents
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    if len(methoxy_matches) < 1:
        return False, "No methoxy groups found"
    
    # Further checks can be done to match additional known rotenoid features (optional)
    # e.g., additional oxygen heterocycles, glycoside linkages, etc.
    
    return True, "Contains chromene-like backbone with methoxy group(s)"

# Example use
# an example SMILES for a rotenoid could be validated here:
# result, reason = is_rotenoid("SMILES HERE")
# print(f"Is rotenoid: {result}, Reason: {reason}")