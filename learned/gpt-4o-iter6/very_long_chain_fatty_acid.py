"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determine if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid has a chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Track visited atoms and longest chain length
    def longest_chain_from(atom_idx, visited):
        max_length = 0
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                length = 1 + longest_chain_from(neighbor.GetIdx(), visited)
                max_length = max(max_length, length)
        visited.remove(atom_idx)
        return max_length

    # Check for the carbon chain attached to an acid, excluding the acid
    chain_lengths = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Only consider carbon atoms
            visited = set()
            chain_length = longest_chain_from(atom.GetIdx(), visited)
            chain_lengths.append(chain_length)

    # Filtering for valid longest chains excluding that which includes the carboxylic acid carbons
    if chain_lengths:
        longest_chain = max(chain_lengths)
    else:
        return False, "No carbon chains found"

    # Determine if the longest carbon chain is a very long-chain fatty acid
    if longest_chain > 22:
        return True, f"Contains {longest_chain} carbons in a continuous chain, qualifies as a very long-chain fatty acid"
    else:
        return False, f"Contains {longest_chain} carbons in the longest chain, does not exceed C22"

# Example usage:
# result, reason = is_very_long_chain_fatty_acid("C(O)(=O)CCCCCCCCCCCCCCCCCCCCCC(O)=O")
# print(result, reason)