"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
from rdkit import Chem

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    A methyl-branched fatty acid is a branched-chain fatty acid containing methyl branches only.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for cyclic structures
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is cyclic"

    # Identify carboxylic acid group (C(=O)O[H])
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)O')
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, f"Found {len(carboxylic_acid_matches)} carboxylic acid groups, need exactly 1"

    # Get the carboxyl carbon atom index
    carboxyl_c_idx = carboxylic_acid_matches[0][0]
    carboxyl_o_indices = [carboxylic_acid_matches[0][1], carboxylic_acid_matches[0][2]]

    # Check for heteroatoms other than oxygens in carboxyl group
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        atom_idx = atom.GetIdx()
        if atomic_num == 6 or atomic_num == 1:
            continue  # Carbon or Hydrogen is allowed
        elif atomic_num == 8:
            if atom_idx in carboxyl_o_indices:
                continue  # Oxygen in carboxyl group is allowed
            else:
                return False, f"Oxygen atom at index {atom_idx} not in carboxylic acid group"
        else:
            return False, f"Heteroatom {atom.GetSymbol()} found at index {atom_idx}"

    # Find the longest carbon chain starting from the carboxyl carbon
    main_chain_atoms = get_longest_chain(mol, carboxyl_c_idx)
    if len(main_chain_atoms) < 2:
        return False, "Main carbon chain is too short"

    # Identify branches and check if they are methyl groups
    for atom_idx in main_chain_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in main_chain_atoms:
                # Check if the branch is a methyl group
                if not is_methyl_branch(mol, neighbor_idx, main_chain_atoms):
                    return False, f"Branch at atom {atom_idx} is larger than methyl group or contains non-carbon atom"

    return True, "Molecule is a methyl-branched fatty acid"

def get_longest_chain(mol, start_idx):
    """
    Finds the longest chain of carbon atoms starting from a given atom index.

    Args:
        mol (Chem.Mol): RDKit molecule object
        start_idx (int): Starting atom index

    Returns:
        list: Atom indices of the longest carbon chain
    """
    max_chain = []
    stack = [(start_idx, [start_idx])]
    visited = set()
    while stack:
        current_idx, path = stack.pop()
        if len(path) > len(max_chain):
            max_chain = path
        atom = mol.GetAtomWithIdx(current_idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if mol.GetAtomWithIdx(neighbor_idx).GetAtomicNum() == 6 and neighbor_idx not in path:
                stack.append((neighbor_idx, path + [neighbor_idx]))
    return max_chain

def is_methyl_branch(mol, branch_start_idx, main_chain_atoms):
    """
    Checks if the branch starting from the given atom index is a methyl group.

    Args:
        mol (Chem.Mol): RDKit molecule object
        branch_start_idx (int): Starting atom index of the branch
        main_chain_atoms (set): Set of atom indices in the main chain

    Returns:
        bool: True if the branch is a methyl group, False otherwise
    """
    visited = set()
    stack = [branch_start_idx]
    branch_atoms = []
    while stack:
        atom_idx = stack.pop()
        if atom_idx in visited or atom_idx in main_chain_atoms:
            continue
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:
            # Branch contains non-carbon atom
            return False
        branch_atoms.append(atom_idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in visited:
                stack.append(neighbor_idx)
    # The branch should contain exactly one carbon atom (methyl group)
    if len(branch_atoms) != 1:
        return False
    return True