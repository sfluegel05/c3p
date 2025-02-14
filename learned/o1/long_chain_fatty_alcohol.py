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
    alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)

    if not alcohol_matches:
        return False, "No alcohol group found"

    # Initialize variables
    min_chain_length = 13
    max_chain_length = 22
    max_found_chain_length = 0

    # For each alcohol group, find the longest carbon chain including the alcohol carbon
    for match in alcohol_matches:
        alcohol_carbon_idx = match[0]
        visited_atoms = set()
        chain_length = longest_chain_length(mol, alcohol_carbon_idx, visited_atoms)
        if chain_length > max_found_chain_length:
            max_found_chain_length = chain_length

    if min_chain_length <= max_found_chain_length <= max_chain_length:
        return True, f"Contains long-chain fatty alcohol with chain length {max_found_chain_length}"
    else:
        return False, f"Longest chain length {max_found_chain_length} does not fall within C{min_chain_length}-C{max_chain_length}"

def longest_chain_length(mol, atom_idx, visited_atoms):
    """
    Recursively explores the molecule to find the longest chain length starting from atom_idx.

    Args:
        mol (Chem.Mol): RDKit molecule
        atom_idx (int): Index of the starting atom
        visited_atoms (set): Set of atom indices that have been visited

    Returns:
        int: Maximum chain length found
    """
    visited_atoms.add(atom_idx)
    atom = mol.GetAtomWithIdx(atom_idx)
    if atom.GetAtomicNum() != 6:
        return 0

    max_length = 1  # Counting the current carbon atom

    for neighbor in atom.GetNeighbors():
        neighbor_idx = neighbor.GetIdx()
        if neighbor_idx not in visited_atoms and neighbor.GetAtomicNum() == 6:
            # Copy visited_atoms to avoid sharing between recursive calls
            length = 1 + longest_chain_length(mol, neighbor_idx, visited_atoms.copy())
            if length > max_length:
                max_length = length

    return max_length