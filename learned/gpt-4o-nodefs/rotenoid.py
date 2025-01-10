"""
Classifies: CHEBI:71543 rotenoid
"""
from rdkit import Chem

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    Rotenoids have distinctive polycyclic systems, often isoflavonoid-like, with specific oxygenation patterns including methoxy groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a rotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define updated SMARTS patterns for rotenoid substructures
    # Isoflavonoid-like or chromene-like backbone with oxygenation
    chromene_pattern = Chem.MolFromSmarts("c1ccc2occc2c1")
    methoxy_pattern = Chem.MolFromSmarts("CO")
    broad_cyclic_structure_pattern = Chem.MolFromSmarts("c1ccc2c(c1)c(=O)c1occc1c2")

    # Check for isoflavonoid or chromene-like structures
    if not mol.HasSubstructMatch(chromene_pattern) and not mol.HasSubstructMatch(broad_cyclic_structure_pattern):
        return False, "Neither chromene-like nor broader isoflavone-like backbone found"

    # Check for methoxy substituents - typically at least one should be present
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    if len(methoxy_matches) < 1:
        return False, "No methoxy groups found"

    # Other potential feature: multiple rings which are part of the backbone
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings() >= 3:
        return False, "Insufficient number of rings for typical rotenoid structure"
    
    return True, "Contains necessary rotenoid structural elements including chromene or isoflavonoid-like backbone with methoxy group(s)"

# Example interaction:
# result, reason = is_rotenoid("SMILES HERE")
# print(f"Is rotenoid: {result}, Reason: {reason}")