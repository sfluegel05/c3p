"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: CHEBI:38839 trienoic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is a polyunsaturated fatty acid with three double bonds,
    typically separated by a methylene group (-CH2-), and a carboxylic acid group
    at the opposite end of the carbon chain from the double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find the carboxylic acid group (-C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"

    carboxyl_atom_idx = carboxyl_matches[0][-1]  # Index of the carbon atom in the carboxyl group

    # Find the double bonds and their positions
    double_bond_atoms = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            double_bond_atoms.extend([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])

    double_bond_atoms = sorted(set(double_bond_atoms))  # Unique atom indices of double-bonded atoms

    # Check if the double bonds are at the opposite end of the carbon chain from the carboxylic acid group
    if not is_carboxyl_opposite_to_double_bonds(mol, carboxyl_atom_idx, double_bond_atoms):
        return False, "Carboxylic acid group not at the opposite end of the carbon chain from the double bonds"

    # Count double bonds
    n_double_bonds = len(double_bond_atoms) // 2  # Each double bond contributes two atom indices
    if n_double_bonds != 3:
        return False, f"Found {n_double_bonds} double bonds, need exactly 3"

    # Check for specific pattern of double bonds separated by -CH2-
    if not has_trienoic_double_bond_pattern(mol, double_bond_atoms):
        return False, "Double bonds not in the expected trienoic pattern"

    # Check for long carbon chain (>= 10 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Carbon chain too short for fatty acid"

    return True, "Contains a carboxylic acid group at the opposite end of the carbon chain from the double bonds, and the double bonds are in the expected trienoic pattern"

def is_carboxyl_opposite_to_double_bonds(mol, carboxyl_atom_idx, double_bond_atoms):
    """
    Checks if the carboxylic acid group (-C(=O)O) is at the opposite end of the carbon chain
    from the double bonds.

    Args:
        mol (Mol): RDKit molecule object
        carboxyl_atom_idx (int): Index of the carbon atom in the carboxyl group
        double_bond_atoms (list): List of atom indices involved in double bonds

    Returns:
        bool: True if the carboxyl group is at the opposite end of the carbon chain from the double bonds,
              False otherwise
    """
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_atom_idx)
    carboxyl_neighbors = [neighbor_atom.GetIdx() for neighbor_atom in carboxyl_atom.GetNeighbors()]

    # Find the longest chain starting from the carboxyl group
    longest_chain = find_longest_chain(mol, carboxyl_atom_idx, exclude_atoms=double_bond_atoms)

    # Check if the longest chain contains any of the double-bonded atoms
    for atom_idx in longest_chain:
        if atom_idx in double_bond_atoms:
            return False

    return True

def find_longest_chain(mol, start_atom_idx, exclude_atoms=None):
    """
    Finds the longest chain of atoms starting from a given atom index, optionally excluding
    specified atoms.

    Args:
        mol (Mol): RDKit molecule object
        start_atom_idx (int): Index of the starting atom
        exclude_atoms (list, optional): List of atom indices to exclude from the chain

    Returns:
        list: List of atom indices in the longest chain
    """
    if exclude_atoms is None:
        exclude_atoms = []

    visited = set()
    longest_chain = []

    def dfs(curr_atom_idx, chain):
        if curr_atom_idx in visited or curr_atom_idx in exclude_atoms:
            return

        visited.add(curr_atom_idx)
        chain.append(curr_atom_idx)

        curr_atom = mol.GetAtomWithIdx(curr_atom_idx)
        neighbors = [neighbor_atom.GetIdx() for neighbor_atom in curr_atom.GetNeighbors()]

        for neighbor_idx in neighbors:
            dfs(neighbor_idx, chain)

        if len(chain) > len(longest_chain):
            longest_chain[:] = chain[:]

        chain.pop()

    dfs(start_atom_idx, [])
    return longest_chain

def has_trienoic_double_bond_pattern(mol, double_bond_atoms):
    """
    Checks if the double bonds in the molecule follow the expected pattern for trienoic fatty acids,
    with double bonds separated by a methylene group (-CH2-).

    Args:
        mol (Mol): RDKit molecule object
        double_bond_atoms (list): List of atom indices involved in double bonds

    Returns:
        bool: True if the double bonds follow the expected trienoic pattern, False otherwise
    """
    # Check if double bonds are separated by -CH2-
    for i in range(len(double_bond_atoms) - 2):
        atom1_idx = double_bond_atoms[i]
        atom2_idx = double_bond_atoms[i + 1]
        atom3_idx = double_bond_atoms[i + 2]

        if not is_methylene_group(mol, atom1_idx, atom2_idx, atom3_idx):
            return False

    return True

def is_methylene_group(mol, atom1_idx, atom2_idx, atom3_idx):
    """
    Checks if the three atoms form a methylene group (-CH2-) separating two double bonds.

    Args:
        mol (Mol): RDKit molecule object
        atom1_idx (int): Index of the first atom
        atom2_idx (int): Index of the second atom
        atom3_idx (int): Index of the third atom

    Returns:
        bool: True if the three atoms form a methylene group, False otherwise
    """
    atom1 = mol.GetAtomWithIdx(atom1_idx)
    atom2 = mol.GetAtomWithIdx(atom2_idx)
    atom3 = mol.GetAtomWithIdx(atom3_idx)

    # Check if atom2 is a carbon with two hydrogens
    if atom2.GetAtomicNum() != 6 or sum(1 for bond in atom2.GetBonds() if bond.GetBondType() == Chem.BondType.SINGLE) != 2:
        return False

    # Check if atom1 and atom3 are double-bonded to atom2
    if not any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in atom2.GetBonds()):
        return False

    return True