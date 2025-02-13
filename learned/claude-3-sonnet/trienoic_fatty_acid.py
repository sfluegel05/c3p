"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: CHEBI:38839 trienoic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is a polyunsaturated fatty acid with three double bonds,
    typically separated by a methylene group (-CH2-), and a carboxylic acid group.

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

    # Find the double bonds and their positions
    double_bond_atoms = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            double_bond_atoms.extend([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])

    double_bond_atoms = sorted(set(double_bond_atoms))  # Unique atom indices of double-bonded atoms

    # Count double bonds
    n_double_bonds = len(double_bond_atoms) // 2  # Each double bond contributes two atom indices
    if n_double_bonds != 3:
        return False, f"Found {n_double_bonds} double bonds, need exactly 3"

    # Check for the specific trienoic double bond pattern
    if not has_trienoic_double_bond_pattern(mol, double_bond_atoms):
        return False, "Double bonds not in the expected trienoic pattern"

    # Check for long carbon chain (>= 10 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Carbon chain too short for fatty acid"

    return True, "Contains a carboxylic acid group and three double bonds in the expected trienoic pattern"

def has_trienoic_double_bond_pattern(mol, double_bond_atoms):
    """
    Checks if the double bonds in the molecule follow the expected pattern for trienoic fatty acids,
    with double bonds separated by a methylene group (-CH2-) or a combination of methylene groups.

    Args:
        mol (Mol): RDKit molecule object
        double_bond_atoms (list): List of atom indices involved in double bonds

    Returns:
        bool: True if the double bonds follow the expected trienoic pattern, False otherwise
    """
    pattern = r'^(C(=C)C([CH2])*C([CH2])*C([CH2])*C(=C)C([CH2])*C([CH2])*C(=C)C([CH2])*C(=O)O)$'
    smiles = Chem.MolToSmiles(mol)

    if Chem.MolFromSmiles(smiles).HasSubstructMatch(Chem.MolFromSmarts(pattern)):
        return True

    return False