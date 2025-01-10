"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.
    An ultra-long-chain fatty acid is a very long-chain fatty acid with a chain length greater than C27.
    Ultra-long-chain fatty acids are monocarboxylic acids with a linear hydrocarbon chain of more than
    27 carbons, possibly including unsaturation and hydroxyl groups, but without significant branching
    or cyclic structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ultra-long-chain fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H1]')
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxy_matches:
        return False, "No carboxylic acid functional group found"

    # Check for exactly one carboxylic acid group
    if len(carboxy_matches) != 1:
        return False, f"Found {len(carboxy_matches)} carboxylic acid groups, expected exactly 1"

    # Build adjacency list of carbon atoms
    adj_list = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            idx = atom.GetIdx()
            adj_list[idx] = []
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:
                    adj_list[idx].append(neighbor.GetIdx())

    # Function to perform depth-first search to find the longest carbon chain
    def dfs(node, visited):
        visited.add(node)
        max_length = 1  # Start with current carbon atom
        for neighbor in adj_list.get(node, []):
            if neighbor not in visited:
                length = 1 + dfs(neighbor, visited.copy())
                if length > max_length:
                    max_length = length
        return max_length

    max_chain_length = 0
    for node in adj_list.keys():
        length = dfs(node, set())
        if length > max_chain_length:
            max_chain_length = length

    if max_chain_length > 27:
        return True, f"Contains linear hydrocarbon chain of {max_chain_length} carbons"
    else:
        return False, f"Longest linear hydrocarbon chain is {max_chain_length} carbons, which is not greater than 27"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'ultra-long-chain fatty acid',
        'definition': 'Any very long-chain fatty acid which has a chain length greater than C27.',
        'parents': []
    }
}