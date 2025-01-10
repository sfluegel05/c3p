"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: CHEBI:XXXXX short-chain fatty acid
"""
from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    A short-chain fatty acid is an aliphatic monocarboxylic acid with a chain length less than C6.
    Oxygen-containing substituents (e.g., hydroxyl and keto groups) are allowed, but other heteroatoms are not permitted.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (monocarboxylic acid)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, f"Found {len(carboxylic_acid_matches)} carboxylic acid groups, need exactly 1"

    # Check for rings (should be acyclic)
    ring_info = mol.GetRingInfo()
    if ring_info and ring_info.NumRings() > 0:
        return False, "Contains rings, should be acyclic"

    # Check for aromaticity (should be aliphatic)
    has_aromatic_atoms = any(atom.GetIsAromatic() for atom in mol.GetAtoms())
    if has_aromatic_atoms:
        return False, "Contains aromatic atoms, should be aliphatic"

    # Check that molecule contains only C, H, and O atoms
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ('C', 'H', 'O'):
            return False, f"Contains heteroatom {atom.GetSymbol()}, not permitted"

    # Remove the carboxylic acid group to analyze the carbon chain
    mol_no_acid = Chem.RWMol(mol)
    acid_match = carboxylic_acid_matches[0]
    # Remove the hydroxyl oxygen
    mol_no_acid.RemoveAtom(acid_match[2])
    # Remove the carbonyl oxygen
    mol_no_acid.RemoveAtom(acid_match[1])
    # Remove the carbonyl carbon (disconnecting the carboxylic acid group)
    mol_no_acid.RemoveAtom(acid_match[0])

    # Get the largest carbon chain length
    atom_indices = [atom.GetIdx() for atom in mol_no_acid.GetAtoms() if atom.GetAtomicNum() == 6]
    if not atom_indices:
        return False, "No carbon chain found after removing carboxylic acid group"

    # Find the longest path consisting of carbon atoms
    max_chain_length = 0
    for atom_idx in atom_indices:
        length = get_longest_chain_length(mol_no_acid, atom_idx, visited=set())
        if length > max_chain_length:
            max_chain_length = length

    if max_chain_length >= 6:
        return False, f"Longest carbon chain length is {max_chain_length}, must be less than 6"

    return True, "Is an aliphatic monocarboxylic acid with less than 6 carbons in the longest chain and only C, H, O atoms"

def get_longest_chain_length(mol, atom_idx, visited):
    """
    Recursively finds the longest carbon chain starting from a given atom.

    Args:
        mol: RDKit molecule object
        atom_idx: Index of the starting atom
        visited: Set of visited atom indices

    Returns:
        int: Length of the longest carbon chain from the starting atom
    """
    visited.add(atom_idx)
    atom = mol.GetAtomWithIdx(atom_idx)
    if atom.GetAtomicNum() != 6:
        return 0

    lengths = [1]  # Include current atom
    for neighbor in atom.GetNeighbors():
        nbr_idx = neighbor.GetIdx()
        if nbr_idx not in visited and neighbor.GetAtomicNum() == 6:
            length = 1 + get_longest_chain_length(mol, nbr_idx, visited.copy())
            lengths.append(length)

    return max(lengths)