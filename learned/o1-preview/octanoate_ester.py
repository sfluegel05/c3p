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

    # Define ester functional group pattern
    ester_pattern = Chem.MolFromSmarts('[$(C(=O)[O;!R])]')  # Non-cyclic ester group
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # Traverse acyl chain for each ester to find if any is octanoate
    for match in ester_matches:
        carbonyl_c_idx = match[0]  # Index of carbonyl carbon atom
        visited_atoms = set()
        acyl_chain_length = count_acyl_chain_carbons(mol, carbonyl_c_idx, visited_atoms)

        if acyl_chain_length == 8:
            return True, "Contains octanoate ester group with 8-carbon acyl chain"
    return False, "No octanoate ester groups with 8-carbon acyl chain found"

def count_acyl_chain_carbons(mol, atom_idx, visited_atoms):
    """
    Recursively counts the number of carbons in the acyl chain starting from the carbonyl carbon.

    Args:
        mol (Chem.Mol): RDKit molecule object
        atom_idx (int): Index of the current atom
        visited_atoms (set): Set of visited atom indices to avoid cycles

    Returns:
        int: Total number of carbons in the acyl chain
    """
    atom = mol.GetAtomWithIdx(atom_idx)
    visited_atoms.add(atom_idx)
    carbon_count = 0

    if atom.GetAtomicNum() == 6:  # Carbon atom
        carbon_count += 1
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in visited_atoms:
                # Skip oxygens to stay on the acyl chain side
                if neighbor.GetAtomicNum() != 8:
                    carbon_count += count_acyl_chain_carbons(mol, neighbor_idx, visited_atoms)
    return carbon_count

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