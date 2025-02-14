"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
from rdkit import Chem

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is a straight-chain C18 polyunsaturated fatty acid
    with two C=C double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group at the terminal
    carboxylic_acid_group = Chem.MolFromSmarts("C(=O)O")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid_group)
    if len(carboxylic_matches) < 1:
        return False, "No terminal carboxylic acid group found"
    
    # Check for the correct number of carbon atoms (exactly 18)
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(carbon_atoms)
    if c_count != 18:
        return False, f"Carbon count is {c_count}, but must be exactly 18"

    # Check for the presence of exactly two C=C double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) != 2:
        return False, f"Found {len(double_bond_matches)} C=C double bonds, need exactly 2"

    # Ensure the chain is straight and does not have any major branchings
    branch_points = sum(1 for atom in carbon_atoms if atom.GetDegree() > 2)
    if branch_points > 0:
        return False, f"Presence of branching in the molecule"

    return True, "The molecule is an octadecadienoic acid with a straight chain, C18, and two C=C double bonds"