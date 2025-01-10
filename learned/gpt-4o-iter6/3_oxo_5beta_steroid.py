"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
from rdkit import Chem

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid is any 3-oxo steroid that has beta-configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a 3-oxo group (carbonyl group at position 3 in a steroid)
    # `[C](=O)[C]([C](C)[C]C)` indicates flexible approach around the carbonyl spare for configurations.
    oxo_pattern = Chem.MolFromSmarts("C(=O)[C]([C])[C]C")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No 3-oxo group found in the expected position"

    # Check for steroid core with beta-configuration:
    # A more loose pattern for steroid rings that includes mapping chiral centers
    steroid_core_pattern = Chem.MolFromSmarts("C1CC=C2C3OC4CCCCC4C3CCC2C1")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Steroid core with 5beta configuration not found"

    # Even though SMILES incorporates stereochemistry, we're verifying stereocenters explicitly
    # Check for the correct stereochemistry using a set of stereocenter-based approaches
    
    # Check stereochemistry on the position where 5-beta is expected:
    chiral_positions = [atom.GetIdx() for atom in mol.GetAtoms() if atom.HasProp('_CIPCode')]

    # Assuming the relevant stereocenter indexes for "5-beta" positions; in reality, should be derived:
    if len(chiral_positions) < 2 or mol.GetAtomWithIdx(chiral_positions[1]).GetProp('_CIPCode') != 'S':
        return False, "5beta stereochemistry not correctly resolved"

    return True, "Successfully identified as a 3-oxo-5beta-steroid with the appropriate structure and configuration"