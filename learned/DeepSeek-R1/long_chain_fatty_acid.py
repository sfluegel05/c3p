"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: long-chain fatty acid (C13-C22)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Get carbonyl and oxygen atoms from the match
    carbonyl_carbon = mol.GetAtomWithIdx(matches[0][0])

    # Find the alpha carbon (connected to carbonyl, not part of the carboxylic acid oxygens)
    alpha_carbon = None
    for neighbor in carbonyl_carbon.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in matches[0][1:]:
            alpha_carbon = neighbor
            break
    if alpha_carbon is None:
        return False, "No alpha carbon found adjacent to carboxylic acid"

    # Helper function to find longest chain from current atom, excluding previous atom
    def get_longest_chain(current_atom, prev_atom):
        max_length = 0
        for neighbor in current_atom.GetNeighbors():
            if neighbor == prev_atom:
                continue
            if neighbor.GetAtomicNum() == 6:  # Only follow carbon atoms
                chain_length = get_longest_chain(neighbor, current_atom)
                if chain_length > max_length:
                    max_length = chain_length
        return 1 + max_length

    # Calculate chain length from alpha carbon, excluding carbonyl
    chain_length = get_longest_chain(alpha_carbon, carbonyl_carbon)
    total_chain_length = chain_length + 1  # Include carbonyl carbon

    # Check chain length
    if 13 <= total_chain_length <= 22:
        return True, f"Long-chain fatty acid with {total_chain_length} carbons"
    else:
        return False, f"Chain length {total_chain_length} is outside C13-C22 range"