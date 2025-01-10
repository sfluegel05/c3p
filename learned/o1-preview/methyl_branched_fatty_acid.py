"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: methyl-branched fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Check for carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Exclude molecules with rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring structures"

    # Find all carboxylic acid groups
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"

    # Assume the molecule is a fatty acid if it contains a carboxylic acid group
    # Identify the carbon atom of the carboxylic acid group
    carboxyl_carbons = [match[0] for match in carboxylic_acid_matches]  # First atom in match is carbon
    # For simplicity, consider the first carboxylic acid group
    carboxyl_c_idx = carboxyl_carbons[0]

    # Find the longest carbon chain ending at the carboxyl carbon
    def get_longest_chain(atom_idx, visited):
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom_idx in visited or atom.GetAtomicNum() != 6:
            return []
        visited.add(atom_idx)
        max_path = []
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if neighbor.GetAtomicNum() == 6 and n_idx != carboxyl_c_idx:
                path = get_longest_chain(n_idx, visited.copy())
                if len(path) > len(max_path):
                    max_path = path
        return [atom_idx] + max_path

    longest_chain = get_longest_chain(carboxyl_c_idx, set())
    if len(longest_chain) < 2:
        return False, "Main carbon chain is too short"

    # Collect atoms in the main chain
    main_chain_atoms = set(longest_chain)

    # Check branches off the main chain
    for atom_idx in longest_chain:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:
            continue
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx not in main_chain_atoms:
                branch_atom = mol.GetAtomWithIdx(n_idx)
                if branch_atom.GetAtomicNum() != 6:
                    return False, "Branch contains heteroatoms"
                # Check if the branch is a methyl group
                branch_size = 1
                branch_visited = {atom_idx}
                atoms_to_visit = [n_idx]
                while atoms_to_visit:
                    current_atom_idx = atoms_to_visit.pop()
                    branch_visited.add(current_atom_idx)
                    current_atom = mol.GetAtomWithIdx(current_atom_idx)
                    for nbr in current_atom.GetNeighbors():
                        nbr_idx = nbr.GetIdx()
                        if nbr_idx not in branch_visited and nbr_idx not in main_chain_atoms:
                            if nbr.GetAtomicNum() != 1:  # Exclude hydrogens
                                branch_size += 1
                                atoms_to_visit.append(nbr_idx)
                if branch_size > 1:
                    return False, f"Branch larger than methyl group found at atom index {atom_idx}"

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