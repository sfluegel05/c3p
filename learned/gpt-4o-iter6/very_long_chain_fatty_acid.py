"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

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
    
    # Check for carboxylic acid group presence
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Use RDKit's FindMolChiralCenters to find the longest chain of continuous carbons
    longest_chain_length = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            visited = set()
            stack = [(atom.GetIdx(), 0)]  # Stack for DFS: (current atom, current chain length)

            while stack:
                current_idx, chain_length = stack.pop()

                if current_idx in visited:
                    continue

                visited.add(current_idx)
                current_atom = mol.GetAtomWithIdx(current_idx)

                # We only count carbon atoms
                if current_atom.GetAtomicNum() == 6:
                    chain_length += 1
                    longest_chain_length = max(longest_chain_length, chain_length)
                    
                # Explore neighbors
                for neighbor in current_atom.GetNeighbors():
                    if neighbor.GetIdx() not in visited:
                        stack.append((neighbor.GetIdx(), chain_length))

    # Determine if the longest carbon chain qualifies
    if longest_chain_length > 22:
        return True, f"Contains {longest_chain_length} carbons in a continuous chain, qualifies as a very long-chain fatty acid"
    else:
        return False, f"Contains {longest_chain_length} carbons in the longest chain, does not exceed C22"