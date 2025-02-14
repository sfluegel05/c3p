"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: CHEBI:36689 trienoic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is a polyunsaturated fatty acid that contains three double bonds.

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

    # Count double bonds
    num_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if num_double_bonds != 3:
        return False, f"Found {num_double_bonds} double bonds, need exactly 3"

    # Check for fatty acid pattern (long carbon chain with carboxylic acid group)
    fatty_acid_pattern = Chem.MolFromSmarts("CCC(=O)O")
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 12:
        return False, "Carbon chain too short for fatty acid"

    # Check for polyunsaturated pattern (double bonds separated by at least one carbon)
    unsaturated_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")
    if not mol.HasSubstructMatch(unsaturated_pattern):
        return False, "Double bonds not separated by at least one carbon"

    # Check for cyclic structures (fatty acids should be acyclic)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C1CCCCC1")):
        return False, "Contains cyclic structures, fatty acids should be acyclic"

    return True, "Contains three double bonds and a carboxylic acid group in a long carbon chain"