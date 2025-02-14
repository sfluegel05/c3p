"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: CHEBI:15904 long-chain fatty acid
"""

from rdkit import Chem

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid based on its SMILES string.
    A long-chain fatty acid is defined as a fatty acid with a chain length ranging from C13 to C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for carboxylic acid group
    carboxylic_acid_smarts = '[CX3](=O)[OX1H0-,OX2H1]'
    carboxylic_acid = Chem.MolFromSmarts(carboxylic_acid_smarts)
    matches = mol.GetSubstructMatches(carboxylic_acid)

    # Check for exactly one carboxylic acid group
    if len(matches) != 1:
        return False, f"Found {len(matches)} carboxylic acid groups, expected exactly 1"

    # Get the carboxyl carbon index and oxygen indices
    carboxyl_match = matches[0]
    carboxyl_carbon_idx = carboxyl_match[0]
    carboxyl_oxygen_idxs = set(carboxyl_match[1:])

    # Exclude carboxyl oxygens from chain search
    excluded_atoms = carboxyl_oxygen_idxs

    # Compute the longest carbon chain starting from the carboxyl carbon
    chain_length = get_longest_chain_length(mol, carboxyl_carbon_idx, excluded_atoms)

    # Check if the chain length is within the specified range
    if 13 <= chain_length <= 22:
        return True, f"Chain length is {chain_length}, within 13-22"
    else:
        return False, f"Chain length is {chain_length}, not within 13-22"

def get_longest_chain_length(mol, start_atom_idx, excluded_atoms):
    """
    Finds the length of the longest carbon chain starting from the start atom.

    Args:
        mol (Mol): RDKit molecule object
        start_atom_idx (int): Atom index to start the search from
        excluded_atoms (set): Set of atom indices to exclude from the search

    Returns:
        int: Length of the longest carbon chain
    """
    max_length = 0
    stack = [(start_atom_idx, 1, {start_atom_idx})]  # (current atom idx, length, visited set)

    while stack:
        current_atom_idx, length, visited = stack.pop()
        if length > max_length:
            max_length = length

        current_atom = mol.GetAtomWithIdx(current_atom_idx)

        for neighbor in current_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()

            if neighbor_idx in visited:
                continue
            if neighbor_idx in excluded_atoms:
                continue
            if mol.GetAtomWithIdx(neighbor_idx).GetAtomicNum() != 6:
                continue  # Only consider carbon atoms

            new_visited = visited.copy()
            new_visited.add(neighbor_idx)
            stack.append((neighbor_idx, length + 1, new_visited))

    return max_length