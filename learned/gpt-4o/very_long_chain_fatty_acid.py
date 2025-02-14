"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid has a chain length greater than 22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure the molecule has a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    carboxyl_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"

    def dfs_chain_length(atom_index, visited, chain_length):
        visited.add(atom_index)
        max_length = chain_length

        for neighbor in mol.GetAtomWithIdx(atom_index).GetNeighbors():
            neighbor_index = neighbor.GetIdx()
            if neighbor_index not in visited and neighbor.GetAtomicNum() == 6:
                max_length = max(max_length, dfs_chain_length(neighbor_index, visited.copy(), chain_length + 1))

        return max_length

    max_chain_length = 0

    # Measure max chain length from each carboxylic acid group's primary carbon atom
    for match in carboxyl_matches:
        carbon_in_carboxyl = match[0]  # The first atom in the match is the carbon atom
        chain_length = dfs_chain_length(carbon_in_carboxyl, set(), 1)  # Start the chain length count from the first carbon

        if chain_length > max_chain_length:
            max_chain_length = chain_length

    # Check chain length criteria for very long-chain fatty acid (greater than 22 carbons)
    if max_chain_length > 22:
        return True, f"Longest carbon chain length is {max_chain_length}, which is greater than 22"

    return False, f"Longest carbon chain length is {max_chain_length}, not greater than 22"