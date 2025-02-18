"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: long-chain fatty acid (C13-C22)
"""
from rdkit import Chem

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

    # Define carboxylic acid pattern (protonated or deprotonated)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1,O-]")
    matches = mol.GetSubstructMatches(carboxylic_acid_pattern)

    # Check for exactly one carboxylic acid group
    if len(matches) != 1:
        return False, f"Found {len(matches)} carboxylic acid groups, need exactly 1"

    # Get carbonyl carbon (first atom in the match)
    carbonyl_carbon = mol.GetAtomWithIdx(matches[0][0])

    # Iterative DFS to find longest carbon chain from carbonyl carbon
    max_length = 1  # Start with carbonyl carbon itself
    stack = [(carbonyl_carbon, None, 1)]  # (current_atom, previous_atom, current_length)
    
    while stack:
        current_atom, prev_atom, length = stack.pop()
        if length > max_length:
            max_length = length
        
        # Explore all neighboring carbons except previous atom
        for neighbor in current_atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 6 or neighbor == prev_atom:
                continue
            
            # Check if the bond is part of a ring to avoid counting ring atoms multiple times
            bond = mol.GetBondBetweenAtoms(current_atom.GetIdx(), neighbor.GetIdx())
            if bond.IsInRing():
                continue  # Skip ring bonds to prevent cyclic path inflation
            
            stack.append((neighbor, current_atom, length + 1))

    # Check chain length
    if 13 <= max_length <= 22:
        return True, f"Long-chain fatty acid with {max_length} carbons"
    elif max_length < 13:
        return False, f"Chain too short ({max_length} carbons)"
    else:
        return False, f"Chain too long ({max_length} carbons)"