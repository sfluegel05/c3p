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

    # Identify alcohol groups (hydroxyl group attached to a carbon)
    alcohol_pattern = Chem.MolFromSmarts("[OX2H][#6]")  # Matches -O(H)-C
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)

    if not alcohol_matches:
        return False, "No alcohol group found"

    min_chain_length = 13
    max_chain_length = 22

    for match in alcohol_matches:
        oxygen_idx = match[0]
        carbon_idx = match[1]

        # Calculate the longest carbon chain starting from the carbon attached to the OH group
        chain_length = get_longest_chain_length(mol, carbon_idx)
        if min_chain_length <= chain_length <= max_chain_length:
            return True, f"Contains long-chain fatty alcohol with chain length {chain_length}"

    return False, "No carbon chain of appropriate length found connected to the alcohol group"


def get_longest_chain_length(mol, start_atom_idx):
    """
    Finds the length of the longest acyclic carbon chain starting from the given atom index.

    Args:
        mol (Chem.Mol): RDKit molecule
        start_atom_idx (int): Index of the starting atom (carbon attached to OH group)

    Returns:
        int: Length of the longest carbon chain (number of carbon atoms)
    """
    max_length = [0]

    def dfs(current_atom_idx, length, visited):
        atom = mol.GetAtomWithIdx(current_atom_idx)
        # Only consider carbon atoms
        if atom.GetAtomicNum() != 6:
            return

        # Update max_length
        if length > max_length[0]:
            max_length[0] = length

        visited.add(current_atom_idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            neighbor_atom = mol.GetAtomWithIdx(neighbor_idx)
            # Continue traversal if neighbor is carbon and not visited
            if neighbor_atom.GetAtomicNum() == 6 and neighbor_idx not in visited:
                dfs(neighbor_idx, length + 1, visited)
        visited.remove(current_atom_idx)

    dfs(start_atom_idx, 1, set())
    return max_length[0]