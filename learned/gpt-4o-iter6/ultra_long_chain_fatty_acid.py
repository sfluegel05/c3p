"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.
    Ultra-long-chain fatty acids have a chain length greater than C27 and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ultra-long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    def dfs_longest_chain(atom_idx, depth, visited):
        """Helper function using DFS to find the longest carbon chain."""
        visited.add(atom_idx)
        max_depth = depth
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                max_depth = max(max_depth, dfs_longest_chain(neighbor.GetIdx(), depth + 1, visited))
        visited.remove(atom_idx)
        return max_depth
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None, "Invalid SMILES string")

    # Ensure the molecule has a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    
    if not matches:
        return False, "No carboxylic acid group found"
    
    # Assuming the carbon of the C(=O) group is the starting point
    carboxyl_carbon_idx = matches[0][0]
    
    # Perform DFS to find the longest carbon chain
    longest_chain = 0
    for neighbor in mol.GetAtomWithIdx(carboxyl_carbon_idx).GetNeighbors():
        if neighbor.GetAtomicNum() == 6:
            longest_chain = max(longest_chain, dfs_longest_chain(neighbor.GetIdx(), 1, set()))

    # Criteria for classification as an ultra-long-chain fatty acid
    if longest_chain > 27:
        return True, f"Contains carboxylic acid group and chain length is C{longest_chain}, which is greater than C27"
    else:
        return False, f"Chain length is C{longest_chain}, which is not greater than C27"