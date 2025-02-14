"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit import RDLogger

# Disable RDKit warnings
RDLogger.DisableLog('rdApp.*')

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    A long-chain fatty alcohol typically has a primary alcohol group and a carbon chain length of 13 to 22 carbons.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify all primary alcohol groups
    primary_alcohol_pattern = Chem.MolFromSmarts("[CX3][OX2H]")
    if not mol.HasSubstructMatch(primary_alcohol_pattern):
        return False, "No primary alcohol group found"
    
    # Get the indices of all alcohol groups containing oxygen
    alcohol_matches = mol.GetSubstructMatches(primary_alcohol_pattern)
    alcohol_indices = {match[1] for match in alcohol_matches} # Collecting oxygen indices
    
    # Use BFS to determine the longest carbon chain from each alcohol group
    max_chain_len = 0
    for alcohol_idx in alcohol_indices:
        visited = set()
        queue = [(alcohol_idx, 0)]
        while queue:
            curr_idx, length = queue.pop(0)
            if curr_idx in visited:
                continue
            visited.add(curr_idx)
            
            atom = mol.GetAtomWithIdx(curr_idx)
            if atom.GetAtomicNum() == 6:  # Count carbon atoms only
                length += 1
            
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited:
                    queue.append((neighbor_idx, length))
        
        # Updating the maximum chain length found
        if max_chain_len < length:
            max_chain_len = length
    
    # Determine if any carbon chain is within the specified length range
    if 13 <= max_chain_len <= 22:
        return True, f"Contains a primary alcohol group and a carbon chain length of {max_chain_len}."
    else:
        return False, f"Longest carbon chain length is {max_chain_len}, outside the 13-22 range."