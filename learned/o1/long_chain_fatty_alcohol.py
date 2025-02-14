"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: Long-chain fatty alcohol
Definition: A fatty alcohol with a chain length ranging from C13 to C22.
"""

from rdkit import Chem

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    A long-chain fatty alcohol has a chain length of 13 to 22 carbons and contains at least one alcohol group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify alcohol groups (hydroxyl groups attached to carbon atoms)
    # Include primary and secondary alcohols
    alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)

    if not alcohol_matches:
        return False, "No alcohol group found"

    min_chain_length = 13
    max_chain_length = 22
    max_found_chain_length = 0

    # For each alcohol group, find the longest carbon chain starting from it
    for match in alcohol_matches:
        alcohol_carbon_idx = match[0]
        chain_length = longest_carbon_chain(mol, alcohol_carbon_idx)
        if chain_length > max_found_chain_length:
            max_found_chain_length = chain_length

    if min_chain_length <= max_found_chain_length <= max_chain_length:
        return True, f"Contains long-chain fatty alcohol with chain length {max_found_chain_length}"
    else:
        return False, f"Longest chain length {max_found_chain_length} does not fall within C{min_chain_length}-C{max_chain_length}"

def longest_carbon_chain(mol, start_atom_idx):
    """
    Finds the length of the longest carbon chain starting from the given atom index.

    Args:
        mol (Chem.Mol): RDKit molecule
        start_atom_idx (int): Index of the starting atom (carbon attached to OH group)

    Returns:
        int: Length of the longest carbon chain (number of carbon atoms)
    """
    visited = set()
    max_length = [0]  # Use a list to allow modification within nested function

    def dfs(atom_idx, length):
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:
            visited.remove(atom_idx)
            return
        # Update max_length if a longer chain is found
        if length > max_length[0]:
            max_length[0] = length
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in visited and mol.GetAtomWithIdx(neighbor_idx).GetAtomicNum() == 6:
                dfs(neighbor_idx, length + 1)
        visited.remove(atom_idx)

    dfs(start_atom_idx, 1)
    return max_length[0]