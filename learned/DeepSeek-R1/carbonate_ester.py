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
    
    # Define SMARTS pattern for the carbonate group O-C(=O)-O
    carbonate_pattern = Chem.MolFromSmarts('[O][C](=O)[O]')
    matches = mol.GetSubstructMatches(carbonate_pattern)
    if not matches:
        return False, "No carbonate (O-C(=O)-O) group detected"
    
    # Check each match for valid carbonate ester structure
    valid = False
    for match in matches:
        # Each match must have exactly 3 atoms (O, C, O)
        if len(match) != 3:
            continue  # Skip invalid matches
        o1_idx, c_idx, o2_idx = match
        
        # Check first oxygen (O1)
        o1 = mol.GetAtomWithIdx(o1_idx)
        if o1.GetTotalNumHs() > 0:
            continue  # Hydrogen present, not ester
        # Check if O1 is connected to a carbon (other than central C)
        o1_neighbors = [n for n in o1.GetNeighbors() if n.GetIdx() != c_idx]
        if not any(n.GetAtomicNum() == 6 for n in o1_neighbors):
            continue  # No organic substituent
        
        # Check second oxygen (O2)
        o2 = mol.GetAtomWithIdx(o2_idx)
        if o2.GetTotalNumHs() > 0:
            continue  # Hydrogen present, not ester
        o2_neighbors = [n for n in o2.GetNeighbors() if n.GetIdx() != c_idx]
        if not any(n.GetAtomicNum() == 6 for n in o2_neighbors):
            continue  # No organic substituent
        
        # If any match passes all checks, return True
        valid = True
        break
    
    if valid:
        return True, "Contains carbonate group with organic substituents"
    else:
        return False, "No valid carbonate ester groups found"