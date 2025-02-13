"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.
    An ultra-long-chain fatty acid is a very long-chain fatty acid with a chain length greater than C27.

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

    # Get the carboxylic acid carbon atom index
    carboxy_c_idx = carboxy_matches[0][0]

    # Initialize variables for DFS
    visited = set()
    max_chain_length = [0]  # Use list for mutable integer in nested function

    def dfs(atom_idx, current_length):
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:
            return
        visited.add(atom_idx)
        # Increment the chain length
        current_length += 1
        # Update max chain length
        if current_length > max_chain_length[0]:
            max_chain_length[0] = current_length
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            bond = mol.GetBondBetweenAtoms(atom_idx, neighbor_idx)
            # Only traverse single or double bonds (to include unsaturations)
            if bond.GetBondType() in (Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE):
                # Avoid cycles
                if neighbor_idx not in visited:
                    dfs(neighbor_idx, current_length)
        visited.remove(atom_idx)

    # Start DFS from the carboxylic acid carbon atom
    dfs(carboxy_c_idx, 0)

    # The chain length includes the carboxylic acid carbon
    chain_length = max_chain_length[0]
    if chain_length > 27:
        return True, f"Longest carbon chain is {chain_length} carbons"

    return False, f"Longest carbon chain is {chain_length} carbons, which is not greater than 27"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'ultra-long-chain fatty acid',
        'definition': 'Any very long-chain fatty acid which has a chain length greater than C27.',
        'parents': []
    }
}