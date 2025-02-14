"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
"""
Classifies: ultra-long-chain fatty acid
"""

from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.
    An ultra-long-chain fatty acid is defined as any very long-chain fatty acid which has a chain length greater than C27.

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

    # Look for carboxylic acid functional group [CX3](=O)[O;H1,-1]
    # This pattern matches both protonated and deprotonated carboxylic acids
    carboxylic_acid = Chem.MolFromSmarts("[CX3](=O)[O;H1,-1]")
    matches = mol.GetSubstructMatches(carboxylic_acid)
    if not matches:
        return False, "No carboxylic acid functional group found"

    max_chain_length = 0

    # For each carboxyl group, find the longest carbon chain
    for match in matches:
        carboxyl_carbon_idx = match[0]
        carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)

        # Function to perform DFS and find the longest chain
        def dfs(atom_idx, visited):
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            max_length = 1  # Include current atom
            lengths = []
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                bond = mol.GetBondBetweenAtoms(atom_idx, neighbor_idx)
                # Only consider carbon atoms not visited yet
                if neighbor_idx not in visited and neighbor.GetAtomicNum() == 6:
                    length = 1 + dfs(neighbor_idx, visited)
                    lengths.append(length)
            visited.remove(atom_idx)
            if lengths:
                return max(lengths)
            else:
                return max_length

        visited = set()
        chain_length = dfs(carboxyl_carbon_idx, visited)
        if chain_length > max_chain_length:
            max_chain_length = chain_length

    if max_chain_length > 27:
        return True, f"Longest carbon chain is {max_chain_length} carbons"
    else:
        return False, f"Longest carbon chain is {max_chain_length} carbons, which is not greater than 27"