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

    # Find the longest carbon chain
    atoms = mol.GetAtoms()
    carbon_count = [0] * len(atoms)
    
    for atom in atoms:
        if atom.GetSymbol() == 'C':
            # Perform breadth-first search to determine the longest chain from each carbon
            queue = [(atom.GetIdx(), 0)]
            visited = set()
            while queue:
                current_idx, count = queue.pop(0)
                if current_idx not in visited:
                    visited.add(current_idx)
                    current_atom = mol.GetAtomWithIdx(current_idx)
                    if current_atom.GetSymbol() == 'C':
                        count += 1
                    carbon_count[current_idx] = max(carbon_count[current_idx], count)
                    for neighbor in current_atom.GetNeighbors():
                        if neighbor.GetIdx() not in visited:
                            queue.append((neighbor.GetIdx(), count))
    
    max_carbon_chain_length = max(carbon_count)
    
    # Classification based on chain length
    if max_carbon_chain_length > 27:
        return True, f"Ultra-long-chain fatty acid with chain length {max_carbon_chain_length}"
    elif max_carbon_chain_length > 22:
        return True, f"Very long-chain fatty acid with chain length {max_carbon_chain_length}"
    else:
        return False, f"Chain length too short, only {max_carbon_chain_length} carbons"

# Example usage:
# print(is_very_long_chain_fatty_acid("OC(=O)CCCCCCCCCCCCCCCCCCCCC=C"))  # Example SMILES