"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
from rdkit import Chem

def is_long_chain_fatty_acid(smiles: str):
    """
    Determine if a molecule is a long-chain fatty acid based on its SMILES string.
    A long-chain fatty acid has a chain length ranging from C13 to C22 with a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Validate the presence of a carboxylic acid group (C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    max_chain_length = 0
    
    # Get all carboxylic acid groups
    acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)

    for match in acid_matches:
        c_idx = match[0]  # Index of carbon in carboxylic group

        # Perform a depth-first search (DFS) to count the longest carbon chain
        visited = set()
        def dfs(atom_idx, chain_length):
            nonlocal max_chain_length
            if chain_length > max_chain_length:
                max_chain_length = chain_length
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:  # only count carbon atoms
                    dfs(neighbor.GetIdx(), chain_length + 1)

        # Start DFS from the carbon attached to the carboxylate group
        dfs(c_idx, 1)  # Start with length 1 as the initial carbon is counted

    if max_chain_length < 13 or max_chain_length > 22:
        return False, f"Longest carbon chain length out of range: {max_chain_length} carbons"
    
    return True, "Contains carboxylic acid group and valid carbon chain length for long-chain fatty acid"