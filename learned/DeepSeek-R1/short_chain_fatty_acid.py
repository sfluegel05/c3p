"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: short-chain fatty acid (CHEBI:26666)
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolFromSmarts

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    A short-chain fatty acid is an aliphatic monocarboxylic acid with a chain length < C6,
    and no non-hydrocarbon substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    mol = MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"

    # Check for exactly one carboxylic acid group
    carboxylic_pattern = MolFromSmarts("[CX3](=O)[OX2H1]")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_pattern)
    if len(carboxylic_matches) != 1:
        return False, f"Found {len(carboxylic_matches)} carboxylic acid groups (needs 1)"

    # Check for other heteroatoms (only C, H, O allowed)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in {1, 6, 8}:
            return False, "Contains non-carbon/hydrogen/oxygen atoms"

    # Get the carbonyl carbon (first atom in the match)
    carbonyl_carbon_idx = carboxylic_matches[0][0]

    # Find longest carbon chain starting from carbonyl carbon
    def dfs(atom_idx, visited, current_length):
        visited.add(atom_idx)
        max_length = current_length
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                length = dfs(neighbor.GetIdx(), visited.copy(), current_length + 1)
                if length > max_length:
                    max_length = length
        return max_length

    chain_length = dfs(carbonyl_carbon_idx, set(), 1)
    if chain_length >= 6:
        return False, f"Main chain length {chain_length} >=6"

    return True, "Aliphatic monocarboxylic acid with chain <6 and no non-hydrocarbon groups"