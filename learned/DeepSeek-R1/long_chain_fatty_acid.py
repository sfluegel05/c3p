"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: long-chain fatty acid (C13-C22)
"""
from rdkit import Chem
from collections import deque

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid (C13-C22) based on its SMILES string.
    A long-chain fatty acid has a carboxylic acid group and an aliphatic chain of 13-22 carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly one carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1,O-]")
    matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(matches) != 1:
        return False, f"Found {len(matches)} carboxylic acid groups, need exactly 1"

    # Get carbonyl carbon and find alpha carbon (adjacent carbon)
    carbonyl_carbon_idx = matches[0][0]
    carbonyl_carbon = mol.GetAtomWithIdx(carbonyl_carbon_idx)
    alpha_carbon = None
    for neighbor in carbonyl_carbon.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:
            alpha_carbon = neighbor
            break
    if not alpha_carbon:
        return False, "No alpha carbon adjacent to carboxylic acid group"

    # BFS to find longest chain starting from alpha_carbon, avoiding revisiting atoms
    max_chain_length = 0
    queue = deque()
    # Each entry is (current_atom_idx, previous_atom_idx, current_length)
    queue.append((alpha_carbon.GetIdx(), carbonyl_carbon_idx, 1))

    while queue:
        current_idx, prev_idx, length = queue.popleft()
        if length > max_chain_length:
            max_chain_length = length
        current_atom = mol.GetAtomWithIdx(current_idx)
        for neighbor in current_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor.GetAtomicNum() != 6 or neighbor_idx == prev_idx:
                continue
            # Proceed only if the neighbor is a carbon and not the previous atom
            queue.append((neighbor_idx, current_idx, length + 1))

    # Total chain length includes the carbonyl carbon (1) + max_chain_length
    total_chain_length = max_chain_length + 1
    if 13 <= total_chain_length <= 22:
        return True, f"Long-chain fatty acid with {total_chain_length} carbons"
    elif total_chain_length < 13:
        return False, f"Chain too short ({total_chain_length} carbons)"
    else:
        return False, f"Chain too long ({total_chain_length} carbons)"