"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid is defined as having a carbon chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Verify carboxylic acid group (-COOH)
    carboxylic_acid = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Correct chain length detection by depth-first search for longest path
    def get_longest_chain_length(mol, start_atom_idx, visited):
        max_length = 0
        stack = [(start_atom_idx, 0)]
        while stack:
            current_idx, length = stack.pop()
            if current_idx not in visited:
                visited.add(current_idx)
                current_atom = mol.GetAtomWithIdx(current_idx)
                if current_atom.GetSymbol() == 'C':
                    length += 1
                for neighbor in current_atom.GetNeighbors():
                    if neighbor.GetIdx() not in visited:
                        stack.append((neighbor.GetIdx(), length))
                max_length = max(max_length, length)
        return max_length

    # Find the longest carbon chain
    max_carbon_chain_length = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            visited = set()
            max_carbon_chain_length = max(max_carbon_chain_length, get_longest_chain_length(mol, atom.GetIdx(), visited))

    # Classification based on chain length with better boundary handling
    if max_carbon_chain_length > 27:
        return True, f"Ultra-long-chain fatty acid with chain length {max_carbon_chain_length}"
    elif max_carbon_chain_length > 22:
        return True, f"Very long-chain fatty acid with chain length {max_carbon_chain_length}"
    elif max_carbon_chain_length == 22:
        # Additional checks can be added for borderline cases, e.g., functional groups
        return True, "Borderline very long-chain fatty acid with exactly 22 carbons, check additional characteristics"
    else:
        return False, f"Chain length too short, only {max_carbon_chain_length} carbons"

# Example usage:
# print(is_very_long_chain_fatty_acid("OC(=O)CCCCCCCCCCCCCCCCCCCCC=C"))  # Example SMILES