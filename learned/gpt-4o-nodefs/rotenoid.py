"""
Classifies: CHEBI:71543 rotenoid
"""
from rdkit import Chem

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    Rotenoids have distinctive polycyclic systems, often isoflavonoid-like, 
    with specific oxygenation patterns including methoxy groups.

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

    # Define more specific SMARTS patterns for rotenoid substructures
    # Rotenoid backbone with oxygenation patterns
    specific_rotenoid_pattern = Chem.MolFromSmarts("c1(c2c3c(ccc4c3c(ccc4c2cc1)O)O)C(=O)O")  
    methoxy_pattern = Chem.MolFromSmarts("CO")

    # Check for the specific rotenoid backbone structure
    if not mol.HasSubstructMatch(specific_rotenoid_pattern):
        return False, "Specific rotenoid backbone structure not found"

    # Check for methoxy substituents - typically at least two should be present in characteristic positions
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    if len(methoxy_matches) < 2:
        return False, f"Insufficient methoxy groups: found {len(methoxy_matches)}"

    # Check for sufficient cyclic structures (at least three rings)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Insufficient number of rings for typical rotenoid structure"
    
    return True, "Contains specific rotenoid structural elements with multiple methoxy groups"

# Example interaction:
# result, reason = is_rotenoid("SMILES HERE")
# print(f"Is rotenoid: {result}, Reason: {reason}")