"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    A short-chain fatty acid is typically a carboxylic acid with an aliphatic chain of 5 or fewer carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid functional group found"

    # Identify the carbon atoms in the carboxylic acid
    matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    longest_chain_length = 0

    for match in matches:
        # Assume the first carbon is the carboxylic acid carbon
        carboxyl_carbon = match[0]
        
        # Perform a breadth-first search from the carboxyl carbon to find the longest aliphatic chain
        visited = set()
        queue = [(carboxyl_carbon, 0)]  # Start from the carboxyl carbon, initial length 0

        while queue:
            current_atom, length = queue.pop(0)

            # Skip if this atom has been visited
            if current_atom in visited:
                continue

            # Mark this atom as visited
            visited.add(current_atom)
            
            # Update the longest chain length
            longest_chain_length = max(longest_chain_length, length)

            # Add neighbors to the queue if they are carbon atoms and not part of the carboxyl group
            for neighbor in mol.GetAtomWithIdx(current_atom).GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited and neighbor.GetIdx() not in match:
                    queue.append((neighbor.GetIdx(), length + 1))

    # A short-chain fatty acid should have the longest aliphatic chain length of 5 or fewer
    if longest_chain_length > 5:
        return False, f"Aliphatic chain too long ({longest_chain_length}), must be 5 or fewer excluding carboxyl carbon"

    return True, "Contains carboxylic acid group with a short aliphatic chain"