"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.
    An ultra-long-chain fatty acid is characterized by a carbon chain length greater than C27.

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

    # Look for carboxylic acid group (COOH)
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Find the longest carbon chain terminating at the carboxylic acid
    def get_longest_chain_length(atom, visited):
        if atom.GetAtomicNum() != 6:
            return 0
        visited.add(atom.GetIdx())
        max_length = 0
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in visited:
                chain_length = get_longest_chain_length(neighbor, visited.copy())
                max_length = max(max_length, chain_length)
        return max_length + 1

    chain_lengths = []
    for match in mol.GetSubstructMatches(carboxylic_pattern):
        visited = set([match[1], match[0]])  # Start from carbon in the C(=O)O group
        if mol.GetAtomWithIdx(match[0]).GetAtomicNum() == 6:  # Ensure starting from carbon
            chain_length = get_longest_chain_length(mol.GetAtomWithIdx(match[0]), visited)
            chain_lengths.append(chain_length)
    
    if not chain_lengths or max(chain_lengths) <= 27:
        return False, f"Longest carbon chain is {max(chain_lengths)}; not greater than 27"

    return True, f"Contains a carbon chain of length {max(chain_lengths)}, qualifying as ultra-long-chain"