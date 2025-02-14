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

    # Identify carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O,H]')
    matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(matches) != 1:
        return False, f"Found {len(matches)} carboxylic acid groups, need exactly 1"
    carboxyl_c_idx = matches[0][0]  # Index of the carboxyl carbon

    # Build main chain starting from carboxyl carbon
    main_chain = []
    visited = set()
    success = build_main_chain(mol, carboxyl_c_idx, visited, main_chain)
    if not success:
        return False, "Failed to build main chain"

    # Check for branches
    for atom_idx in main_chain:
        atom = mol.GetAtomWithIdx(atom_idx)
        neighbors = atom.GetNeighbors()
        for neighbor in neighbors:
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in main_chain:
                # This is a branch
                branch_size = count_carbon_atoms(neighbor, set([atom_idx]))
                if branch_size != 1:
                    return False, f"Branch at atom {atom_idx + 1} is larger than a methyl group"

    # Check for heteroatoms other than carboxyl group
    carboxyl_oxygen_indices = [idx for match in matches for idx in match if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 6 or atomic_num == 1:
            continue  # Carbon or Hydrogen is allowed
        elif atomic_num == 8:
            # Oxygen atoms should only be in the carboxyl group
            if atom.GetIdx() not in carboxyl_oxygen_indices:
                return False, f"Extra oxygen atom found at index {atom.GetIdx() + 1}"
        else:
            return False, f"Heteroatom {atom.GetSymbol()} found at index {atom.GetIdx() + 1}"

    return True, "Molecule is a methyl-branched fatty acid"


def build_main_chain(mol, atom_idx, visited, chain):
    """
    Recursively builds the main carbon chain starting from the carboxyl carbon.

    Args:
        mol (Chem.Mol): RDKit molecule object
        atom_idx (int): Current atom index
        visited (set): Set of visited atom indices
        chain (list): List to store the main chain atom indices

    Returns:
        bool: True if a terminal carbon is reached, False otherwise
    """
    atom = mol.GetAtomWithIdx(atom_idx)
    visited.add(atom_idx)
    chain.append(atom_idx)
    neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited]

    if len(neighbors) == 0:
        # Reached terminal carbon
        return True
    elif len(neighbors) == 1:
        # Continue along the chain
        return build_main_chain(mol, neighbors[0].GetIdx(), visited, chain)
    else:
        # Multiple paths, choose the one leading to the longest chain
        longest_chain = []
        for nbr in neighbors:
            temp_chain = []
            temp_visited = visited.copy()
            temp_success = build_main_chain(mol, nbr.GetIdx(), temp_visited, temp_chain)
            if temp_success and len(temp_chain) > len(longest_chain):
                longest_chain = temp_chain
        if longest_chain:
            chain.extend(longest_chain)
            return True
    return False


def count_carbon_atoms(atom, visited):
    """
    Counts the number of carbon atoms in a branch recursively.

    Args:
        atom (Chem.Atom): Starting atom
        visited (set): Set of visited atom indices

    Returns:
        int: Number of carbon atoms in the branch
    """
    count = 0
    if atom.GetAtomicNum() == 6:
        count += 1
    else:
        return 0  # Non-carbon atom found in branch
    visited.add(atom.GetIdx())
    for neighbor in atom.GetNeighbors():
        if neighbor.GetIdx() not in visited:
            count += count_carbon_atoms(neighbor, visited)
    return count