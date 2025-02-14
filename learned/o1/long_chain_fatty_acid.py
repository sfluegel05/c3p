"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: CHEBI:15904 long-chain fatty acid
"""

from rdkit import Chem

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid based on its SMILES string.
    A long-chain fatty acid is defined as a fatty acid with a chain length ranging from C13 to C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find carboxylic acid group using SMARTS pattern
    carboxylic_acid_smarts = '[CX3](=O)[OX1H]'
    carboxylic_acid = Chem.MolFromSmarts(carboxylic_acid_smarts)
    matches = mol.GetSubstructMatches(carboxylic_acid)

    if not matches:
        return False, "No carboxylic acid group found"

    # Check that there is exactly one carboxylic acid group
    if len(matches) != 1:
        return False, f"Found {len(matches)} carboxylic acid groups, expected 1"

    # Get the index of the carboxyl carbon atom
    carboxyl_carbon_idx = matches[0][0]

    # Get the longest carbon chain length starting from the carboxyl carbon
    chain_length = get_longest_chain_length(mol, carboxyl_carbon_idx)

    # Check if the chain length is within the specified range
    if 13 <= chain_length <= 22:
        return True, f"Chain length is {chain_length}, within 13-22"
    else:
        return False, f"Chain length is {chain_length}, not within 13-22"

def get_longest_chain_length(mol, start_atom_idx):
    """
    Finds the length of the longest carbon chain starting from the given atom index.

    Args:
        mol (Mol): RDKit molecule object
        start_atom_idx (int): Atom index to start from

    Returns:
        int: Length of the longest carbon chain
    """
    max_length = [0]  # List to allow modification within nested function

    def dfs(current_atom_idx, visited, length):
        """
        Depth-first search to find the longest path of carbon atoms.

        Args:
            current_atom_idx (int): Current atom index
            visited (set): Set of visited atom indices
            length (int): Current chain length
        """
        atom = mol.GetAtomWithIdx(current_atom_idx)
        visited.add(current_atom_idx)
        neighbors = [
            neighbor for neighbor in atom.GetNeighbors()
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited
        ]

        if not neighbors:
            # If no unvisited carbon neighbors, update max_length if necessary
            if length > max_length[0]:
                max_length[0] = length
        else:
            for neighbor in neighbors:
                dfs(neighbor.GetIdx(), visited, length + 1)
        
        visited.remove(current_atom_idx)

    # Start DFS from the carboxyl carbon
    dfs(start_atom_idx, set(), 1)
    return max_length[0]