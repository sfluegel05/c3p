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
    which has an 8-carbon unbranched, acyclic chain including the carbonyl carbon.

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

    # Define SMARTS pattern for octanoate ester
    octanoate_ester_pattern = Chem.MolFromSmarts('C(=O)O[C;!R]')
    ester_matches = mol.GetSubstructMatches(octanoate_ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # Define the octanoyl group as SMARTS (linear chain of 8 carbons including carbonyl carbon)
    octanoyl_smarts = '[C;X3](=O)[C;X4H2][C;X4H2][C;X4H2][C;X4H2][C;X4H2][C;X4H2][C;X4H3]'
    octanoyl_pattern = Chem.MolFromSmarts(octanoyl_smarts)

    # For each ester group, check if the acyl chain matches the octanoyl pattern
    for match in ester_matches:
        carbonyl_c_idx = match[0]  # Index of carbonyl carbon atom

        # Get the acyl chain starting from the carbonyl carbon
        acyl_chain_atom_idxs = get_acyl_chain_atoms(mol, carbonyl_c_idx)

        if acyl_chain_atom_idxs is None:
            continue  # Skip if acyl chain is branched or cyclic

        # Create a sub-molecule of the acyl chain
        acyl_chain = Chem.PathToSubmol(mol, acyl_chain_atom_idxs)

        # Check if acyl chain matches the octanoyl pattern
        if acyl_chain.HasSubstructMatch(octanoyl_pattern):
            return True, "Contains octanoate ester group with 8-carbon unbranched acyl chain"

    return False, "No octanoate ester groups with unbranched 8-carbon acyl chain found"

def get_acyl_chain_atoms(mol, carbonyl_atom_idx):
    """
    Retrieves the atom indices of the acyl chain starting from the carbonyl carbon,
    proceeding along a linear, unbranched, acyclic path.

    Args:
        mol (Chem.Mol): RDKit molecule object
        carbonyl_atom_idx (int): Index of the carbonyl carbon atom

    Returns:
        list: Atom indices of the acyl chain, or None if not linear or acyclic
    """
    acyl_chain_atom_idxs = [carbonyl_atom_idx]
    visited_atoms = set(acyl_chain_atom_idxs)
    current_atom_idx = carbonyl_atom_idx

    # Find the alpha carbon (carbon attached to carbonyl carbon)
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_atom_idx)
    neighbors = [nbr for nbr in carbonyl_atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and not nbr.IsInRing()]
    if len(neighbors) != 1:
        return None  # Acyl chain is branched or cyclic
    current_atom = neighbors[0]
    current_atom_idx = current_atom.GetIdx()
    acyl_chain_atom_idxs.append(current_atom_idx)
    visited_atoms.add(current_atom_idx)

    # Traverse the acyl chain
    while True:
        neighbors = [nbr for nbr in current_atom.GetNeighbors()
                     if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited_atoms and not nbr.IsInRing()]
        if len(neighbors) == 0:
            # Reached end of chain (terminal carbon)
            break
        elif len(neighbors) > 1:
            # Branching detected
            return None  # Not a linear unbranched chain
        else:
            # Proceed to next carbon
            current_atom = neighbors[0]
            current_atom_idx = current_atom.GetIdx()
            acyl_chain_atom_idxs.append(current_atom_idx)
            visited_atoms.add(current_atom_idx)

    # Check if the acyl chain has exactly 8 carbons (including carbonyl carbon)
    if len(acyl_chain_atom_idxs) != 8:
        return None

    return acyl_chain_atom_idxs

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