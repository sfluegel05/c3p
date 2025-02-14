"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: CHEBI:51766 long-chain fatty alcohol
A fatty alcohol with a chain length ranging from C13 to C22.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    A long-chain fatty alcohol has a chain length ranging from C13 to C22 with an -OH group attached.

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

    # Look for an -OH group
    oh_pattern = Chem.MolFromSmarts("[OX1H]")
    oh_indices = mol.GetSubstructMatches(oh_pattern)
    if len(oh_indices) != 1:
        return False, "Molecule should have exactly one hydroxy (-OH) group"

    # Find the longest aliphatic chain
    longest_chain_length = 0
    longest_chain_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            chain_atoms = _find_aliphatic_chain(mol, atom.GetIdx())
            chain_length = len(chain_atoms)
            if chain_length > longest_chain_length:
                longest_chain_length = chain_length
                longest_chain_atoms = chain_atoms

    if not (13 <= longest_chain_length <= 22):
        return False, f"Longest aliphatic chain length is {longest_chain_length}, not in the range C13 to C22"

    # Check if -OH group is attached to the longest aliphatic chain
    oh_atom_idx = oh_indices[0]
    oh_atom = mol.GetAtomWithIdx(oh_atom_idx)
    if oh_atom_idx not in longest_chain_atoms:
        return False, "Hydroxy (-OH) group not attached to the longest aliphatic chain"

    # Check for disallowed functional groups
    disallowed_patterns = [
        Chem.MolFromSmarts("[C$(C=O)O]"),  # Esters
        Chem.MolFromSmarts("C(=O)C"),  # Ketones
        Chem.MolFromSmarts("C(=O)O"),  # Carboxylic acids
    ]
    for pattern in disallowed_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Contains disallowed functional groups"

    return True, "Molecule contains a straight-chain aliphatic alcohol with a chain length between C13 and C22"

def _find_aliphatic_chain(mol, start_atom_idx):
    """
    Finds the aliphatic chain starting from the given atom index.

    Args:
        mol (Mol): RDKit molecule object
        start_atom_idx (int): Index of the starting atom

    Returns:
        list: List of atom indices in the aliphatic chain
    """
    visited = set()
    queue = [start_atom_idx]
    chain_atoms = []

    while queue:
        atom_idx = queue.pop(0)
        if atom_idx in visited:
            continue
        visited.add(atom_idx)

        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:  # Carbon atom
            chain_atoms.append(atom_idx)
            for neighbor_atom in atom.GetNeighbors():
                neighbor_idx = neighbor_atom.GetIdx()
                if neighbor_idx not in visited:
                    queue.append(neighbor_idx)

    return chain_atoms