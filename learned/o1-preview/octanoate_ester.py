"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: octanoate ester
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester is an ester where the acid component is octanoic acid (caprylic acid),
    which has an 8-carbon unbranched chain including the carbonyl carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester functional group pattern (non-cyclic ester)
    ester_pattern = Chem.MolFromSmarts('[$(C(=O)[O;!R])]')  # Non-cyclic ester group
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # For each ester group, check if the acyl chain is unbranched and 8 carbons long
    for match in ester_matches:
        carbonyl_c_idx = match[0]  # Index of carbonyl carbon atom

        acyl_chain_length = count_acyl_chain_length(mol, carbonyl_c_idx)
        if acyl_chain_length == 8:
            return True, "Contains octanoate ester group with 8-carbon unbranched acyl chain"

    return False, "No octanoate ester groups with unbranched 8-carbon acyl chain found"

def count_acyl_chain_length(mol, carbonyl_atom_idx):
    """
    Counts the number of carbons in the acyl chain starting from the carbonyl carbon,
    proceeding along a linear, unbranched path.

    Args:
        mol (Chem.Mol): RDKit molecule object
        carbonyl_atom_idx (int): Index of the carbonyl carbon atom

    Returns:
        int: Total number of carbons in the acyl chain, or None if not linear or acyclic
    """
    length = 1  # Start with the carbonyl carbon
    prev_atom_idx = carbonyl_atom_idx
    current_atom_idx = None

    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_atom_idx)

    # Find the alpha carbon (carbon attached to carbonyl carbon)
    for neighbor in carbonyl_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # Carbon atom
            current_atom_idx = neighbor.GetIdx()
            current_atom = neighbor
            break
    else:
        # No carbon neighbor found
        return None

    visited_atoms = set()
    visited_atoms.add(prev_atom_idx)
    visited_atoms.add(current_atom_idx)

    while True:
        length += 1
        neighbor_atoms = [nbr for nbr in current_atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev_atom_idx]

        if len(neighbor_atoms) == 0:
            # Reached end of chain (terminal carbon)
            break
        elif len(neighbor_atoms) > 1:
            # Branching detected
            return None  # Not a linear unbranched chain
        else:
            # Proceed to next carbon
            prev_atom_idx = current_atom_idx
            current_atom_idx = neighbor_atoms[0].GetIdx()
            if current_atom_idx in visited_atoms:
                # Cycle detected
                return None
            visited_atoms.add(current_atom_idx)
            current_atom = neighbor_atoms[0]

    return length

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'octanoate ester',
        'definition': 'Any fatty acid ester in which the carboxylic acid component is octanoic acid (caprylic acid).',
        'parents': []
    },
    'config': {},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}