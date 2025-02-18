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
    A carbonate ester is carbonic acid where both hydrogens are replaced by organyl groups.

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
    
    # SMARTS pattern for O-C(=O)-O group
    carbonate_pattern = Chem.MolFromSmarts('[O][C](=O)[O]')
    matches = mol.GetSubstructMatches(carbonate_pattern)
    if not matches:
        return False, "No carbonate (O-C(=O)-O) group detected"
    
    # Check each match for organic substituents on both oxygens
    valid = False
    for match in matches:
        if len(match) != 3:
            continue
        o1_idx, c_idx, o2_idx = match
        
        # Check O1's neighbors (excluding central C)
        o1 = mol.GetAtomWithIdx(o1_idx)
        o1_neighbors = [n for n in o1.GetNeighbors() if n.GetIdx() != c_idx]
        if not any(n.GetAtomicNum() == 6 for n in o1_neighbors):
            continue
        
        # Check O2's neighbors (excluding central C)
        o2 = mol.GetAtomWithIdx(o2_idx)
        o2_neighbors = [n for n in o2.GetNeighbors() if n.GetIdx() != c_idx]
        if not any(n.GetAtomicNum() == 6 for n in o2_neighbors):
            continue
        
        valid = True
        break
    
    if valid:
        return True, "Contains carbonate group with organic substituents on both oxygens"
    else:
        return False, "No valid carbonate ester groups found"