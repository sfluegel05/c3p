"""
Classifies: CHEBI:71543 rotenoid
"""
from rdkit import Chem

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    Rotenoids are characterized by a complex polycyclic structure and
    specific oxygenation patterns.

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

    # Refined SMARTS patterns for rotenoid-like substructures
    # Basic isoflavonoid-like core common to many rotenoids
    isoflavonoid_core_pattern = Chem.MolFromSmarts("c1cc2c(cc1)Oc3c(ccc(c3)O2)C(=O)")

    # Methoxy pattern
    methoxy_pattern = Chem.MolFromSmarts("CO")

    # Check for the refined rotenoid-like backbone structure
    if not mol.HasSubstructMatch(isoflavonoid_core_pattern):
        return False, "Rotenoid-like core structure not found"

    # Check for minimum methoxy group presence
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    if len(methoxy_matches) < 2:
        return False, f"Insufficient methoxy groups: found {len(methoxy_matches)}"

    # Additional checks: At least 3 rings are expected in diketone core for rotenoids
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Insufficient number of rings for typical rotenoid structure"

    return True, "Structural features consistent with rotenoids including methoxy groups"

# Example interaction:
# result, reason = is_rotenoid("SMILES HERE")
# print(f"Is rotenoid: {result}, Reason: {reason}")