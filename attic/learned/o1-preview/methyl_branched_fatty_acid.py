"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: methyl-branched fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    A methyl-branched fatty acid is a fatty acid containing only methyl branches.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group at terminal position
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Check for rings (should be acyclic)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring structures"

    # Identify the carboxyl carbon atom
    matches = mol.GetSubstructMatches(carboxylic_acid)
    carboxyl_carbons = [match[0] for match in matches]  # First atom in pattern is the carbon
    if not carboxyl_carbons:
        return False, "No carboxyl carbon found"
    carboxyl_c = carboxyl_carbons[0]

    # Find the longest path starting from the carboxyl carbon
    def find_longest_path(atom_idx, visited):
        visited = visited + [atom_idx]
        paths = [visited]
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if neighbor.GetAtomicNum() != 6:
                continue  # Skip non-carbon atoms
            if n_idx not in visited:
                new_paths = find_longest_path(n_idx, visited)
                if len(new_paths) > 0:
                    paths.extend(new_paths)
        return paths

    all_paths = find_longest_path(carboxyl_c, [])
    longest_path = max(all_paths, key=len)
    main_chain_atoms = set(longest_path)

    # Check for branches off the main chain
    for atom_idx in main_chain_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:
            continue  # Skip non-carbon atoms
        neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        non_chain_neighbors = [n_idx for n_idx in neighbors if n_idx not in main_chain_atoms]
        # If there are branches
        for n_idx in non_chain_neighbors:
            submol_atoms = Chem.rdmolops.GetShortestPath(mol, atom_idx, n_idx)
            branch_atoms = set(submol_atoms) - main_chain_atoms
            # Check if branch consists of only one carbon (methyl group)
            if len(branch_atoms) != 1:
                return False, "Branch larger than methyl group found"
            branch_atom = mol.GetAtomWithIdx(n_idx)
            # Ensure that the branch atom has no further substituents (except hydrogens)
            if any(neighbor.GetAtomicNum() != 1 for neighbor in branch_atom.GetNeighbors() if neighbor.GetIdx() != atom_idx):
                return False, "Branch is not a methyl group"
    return True, "Molecule is a methyl-branched fatty acid with only methyl branches"

__metadata__ = {
    'chemical_class': {
        'name': 'methyl-branched fatty acid',
        'definition': 'Any branched-chain fatty acid containing methyl branches only.'
    },
    'message': None,
    'success': True,
    'error': '',
}