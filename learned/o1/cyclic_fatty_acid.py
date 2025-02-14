"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
"""
Classifies: cyclic fatty acid
"""

from rdkit import Chem

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    A cyclic fatty acid is any fatty acid containing anywhere in its structure a ring of atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for terminal carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No terminal carboxylic acid group found"

    # Assume only one terminal carboxylic acid group
    carboxy_carbon_idx = carboxylic_acid_matches[0][0]  # Carbonyl carbon index

    # Find the longest aliphatic carbon chain connected to the carboxyl carbon
    max_chain_length = get_longest_aliphatic_chain(mol, carboxy_carbon_idx)
    
    # Check if the chain is long enough to be considered a fatty acid (e.g., at least 8 carbons)
    MIN_CHAIN_LENGTH = 8
    if max_chain_length < MIN_CHAIN_LENGTH:
        return False, f"Aliphatic chain too short ({max_chain_length} carbons), not a fatty acid"

    # Check for ring structures
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return False, "No ring structures found"

    return True, f"Contains {num_rings} ring(s) and a fatty acid chain of {max_chain_length} carbons"

def get_longest_aliphatic_chain(mol, carboxy_carbon_idx):
    """
    Finds the longest aliphatic carbon chain connected to the carboxyl carbon.

    Args:
        mol (Chem.Mol): RDKit molecule object
        carboxy_carbon_idx (int): Index of the carboxyl carbon atom

    Returns:
        int: Length of the longest aliphatic carbon chain
    """
    from collections import deque

    max_length = 0

    # Exclude the carboxyl group atoms
    excluded_atoms = set()
    excluded_atoms.add(carboxy_carbon_idx)
    for neighbor in mol.GetAtomWithIdx(carboxy_carbon_idx).GetNeighbors():
        excluded_atoms.add(neighbor.GetIdx())

    # Start from atoms connected to the carboxyl carbon (excluding carboxyl oxygen)
    for neighbor in mol.GetAtomWithIdx(carboxy_carbon_idx).GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # Carbon atom
            start_idx = neighbor.GetIdx()
            visited = set()
            queue = deque()
            queue.append((start_idx, 1))  # (atom index, chain length)
            visited.add(start_idx)

            while queue:
                current_idx, length = queue.popleft()
                current_atom = mol.GetAtomWithIdx(current_idx)
                
                # Only consider aliphatic carbon atoms
                if current_atom.GetAtomicNum() != 6 or current_atom.IsAromatic():
                    continue

                max_length = max(max_length, length)
                
                for nbr in current_atom.GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx in visited or nbr_idx in excluded_atoms:
                        continue
                    visited.add(nbr_idx)
                    queue.append((nbr_idx, length + 1))
    
    return max_length