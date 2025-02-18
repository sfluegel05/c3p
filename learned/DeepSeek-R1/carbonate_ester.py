"""
Classifies: CHEBI:46722 carbonate ester
"""
"""
Classifies: CHEBI:48349 carbonate ester
"""
from rdkit import Chem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    A carbonate ester has a central carbonyl group connected to two oxygen atoms,
    each of which is bonded to an organyl group (no hydrogen substituents).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbonate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for the carbonate group O-C(=O)-O
    carbonate_pattern = Chem.MolFromSmarts('[O][C](=O)[O]')
    matches = mol.GetSubstructMatches(carbonate_pattern)
    if not matches:
        return False, "No carbonate (O-C(=O)-O) group detected"
    
    # Check each oxygen in the carbonate group for organic substituents
    for match in matches:
        # match indices: [O1, C, O2]
        o1_idx, c_idx, o2_idx = match
        
        # Check first oxygen (O1)
        o1 = mol.GetAtomWithIdx(o1_idx)
        if o1.GetTotalNumHs() > 0:
            return False, "Hydrogen present on carbonate oxygen"
        # Check if O1 is connected to a carbon (other than the central C)
        o1_neighbors = [n for n in o1.GetNeighbors() if n.GetIdx() != c_idx]
        if not any(n.GetAtomicNum() == 6 for n in o1_neighbors):
            return False, "Carbonate oxygen not attached to organic group"
        
        # Check second oxygen (O2)
        o2 = mol.GetAtomWithIdx(o2_idx)
        if o2.GetTotalNumHs() > 0:
            return False, "Hydrogen present on carbonate oxygen"
        o2_neighbors = [n for n in o2.GetNeighbors() if n.GetIdx() != c_idx]
        if not any(n.GetAtomicNum() == 6 for n in o2_neighbors):
            return False, "Carbonate oxygen not attached to organic group"
    
    return True, "Contains carbonate group with organic substituents"