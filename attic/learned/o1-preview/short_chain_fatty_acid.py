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

    # Check that molecule contains only C, H, and O atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6, 8):
            return False, f"Contains disallowed atom {atom.GetSymbol()}, only C, H, O are permitted"

    # Check for carboxylic acid group (monocarboxylic acid)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, f"Found {len(carboxylic_acid_matches)} carboxylic acid groups, need exactly 1"

    # Check for rings (should be acyclic)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains rings, should be acyclic"

    # Check for aromaticity (should be aliphatic)
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Contains aromatic atoms, should be aliphatic"

    # Get the carboxyl carbon index
    carboxyl_carbon_idx = carboxylic_acid_matches[0][0]
    carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)

    # Traverse the molecule to find the longest carbon chain starting from the carboxyl carbon
    visited = set()
    chain_length = get_longest_carbon_chain_length(mol, carboxyl_carbon_idx, visited)

    if chain_length >= 6:
        return False, f"Chain length is {chain_length}, must be less than 6 carbons (including carboxyl carbon)"

    return True, "Is an aliphatic monocarboxylic acid with chain length less than C6 and only C, H, O atoms"

def get_longest_carbon_chain_length(mol, atom_idx, visited):
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

    max_length = 1  # Include current atom
    for neighbor in atom.GetNeighbors():
        nbr_idx = neighbor.GetIdx()
        if nbr_idx not in visited and neighbor.GetAtomicNum() == 6:
            length = 1 + get_longest_carbon_chain_length(mol, nbr_idx, visited.copy())
            if length > max_length:
                max_length = length

    return max_length