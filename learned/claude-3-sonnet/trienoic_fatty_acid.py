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
    at the end of the carbon chain.

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
    
    # Check for carboxylic acid group (-C(=O)O) at the end of the carbon chain
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"
    
    carboxyl_atom_idx = carboxyl_matches[0][-1]  # Index of the carbon atom in the carboxyl group
    if not is_terminal_carboxyl(mol, carboxyl_atom_idx):
        return False, "Carboxylic acid group not at the end of the carbon chain"
    
    # Count double bonds
    n_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if n_double_bonds != 3:
        return False, f"Found {n_double_bonds} double bonds, need exactly 3"
    
    # Check for specific pattern of double bonds separated by -CH2-
    if not has_trienoic_double_bond_pattern(mol):
        return False, "Double bonds not in the expected trienoic pattern"
    
    # Check for long carbon chain (>= 10 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Carbon chain too short for fatty acid"
    
    return True, "Contains a carboxylic acid group at the end of the carbon chain and 3 double bonds in the expected trienoic pattern"

def is_terminal_carboxyl(mol, carboxyl_atom_idx):
    """
    Checks if the carboxylic acid group (-C(=O)O) is at the end of the carbon chain.

    Args:
        mol (Mol): RDKit molecule object
        carboxyl_atom_idx (int): Index of the carbon atom in the carboxyl group

    Returns:
        bool: True if the carboxyl group is at the end of the carbon chain, False otherwise
    """
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_atom_idx)
    for neighbor_atom in carboxyl_atom.GetNeighbors():
        if neighbor_atom.GetAtomicNum() == 6:
            return sum(1 for bond in neighbor_atom.GetBonds() if bond.GetBondType() == Chem.BondType.SINGLE) == 1
    return False

def has_trienoic_double_bond_pattern(mol):
    """
    Checks if the double bonds in the molecule follow the expected pattern for trienoic fatty acids,
    with double bonds separated by a methylene group (-CH2-).

    Args:
        mol (Mol): RDKit molecule object

    Returns:
        bool: True if the double bonds follow the expected trienoic pattern, False otherwise
    """
    double_bond_atoms = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            double_bond_atoms.extend([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
    
    double_bond_atoms = sorted(set(double_bond_atoms))  # Unique atom indices of double-bonded atoms
    
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