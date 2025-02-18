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
    matches = mol.GetSubstructMatches(carboxylic_pattern)
    if len(matches) != 1:
        return False, f"Found {len(matches)} carboxylic acid groups (needs 1)"

    # Check all atoms are only C, H, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in {1, 6, 8}:
            return False, "Contains non-C/H/O atoms"

    # Verify exactly two oxygen atoms (carboxylic group only)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 2:
        return False, f"Found {o_count} oxygen atoms (must be exactly 2)"

    # Get carbonyl carbon index from match
    carbonyl_idx = matches[0][0]
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)

    # Find alpha carbon (must be directly bonded to carbonyl)
    alpha_candidates = []
    for neighbor in carbonyl_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:
            alpha_candidates.append(neighbor.GetIdx())
    
    if not alpha_candidates:
        return False, "No alpha carbon found"

    alpha_idx = alpha_candidates[0]

    # DFS to find longest chain from alpha carbon (excluding carbonyl)
    def dfs(current, visited):
        max_length = 1  # Start counting from current atom
        visited.add(current)
        atom = mol.GetAtomWithIdx(current)
        for neighbor in atom.GetNeighbors():
            nbr_idx = neighbor.GetIdx()
            if nbr_idx == carbonyl_idx:
                continue  # Skip back to carbonyl
            if neighbor.GetAtomicNum() == 6 and nbr_idx not in visited:
                length = dfs(nbr_idx, visited.copy())
                max_length = max(max_length, 1 + length)
        return max_length

    chain_length = dfs(alpha_idx, set())

    if chain_length >= 6:
        return False, f"Main chain length {chain_length} >=6"

    return True, "Aliphatic monocarboxylic acid with chain <6 and no non-hydrocarbon groups"