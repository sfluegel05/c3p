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

    # Look for carboxylic acid functional group [CX3](=O)[OX1H]
    carboxylic_acid = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    matches = mol.GetSubstructMatches(carboxylic_acid)
    if not matches:
        return False, "No carboxylic acid functional group found"

    # Assume the first match is the carboxyl carbon
    carboxyl_carbon_idx = matches[0][0]

    # Function to find the longest carbon chain starting from the carboxyl carbon
    def find_longest_chain(atom_idx, visited):
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        max_length = 1  # Include current atom
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in visited and neighbor.GetAtomicNum() == 6:
                length = 1 + find_longest_chain(neighbor_idx, visited)
                if length > max_length:
                    max_length = length
        visited.remove(atom_idx)
        return max_length

    # Start the search from the carboxyl carbon
    visited = set()
    chain_length = find_longest_chain(carboxyl_carbon_idx, visited)
    if chain_length > 27:
        return True, f"Longest carbon chain is {chain_length} carbons"
    else:
        return False, f"Longest carbon chain is {chain_length} carbons, which is not greater than 27"