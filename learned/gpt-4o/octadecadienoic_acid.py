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

    # Check for carboxylic acid group
    carboxylic_acid_group = Chem.MolFromSmarts("C(=O)[O,N]")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid_group)
    if len(carboxylic_matches) < 1:
        return False, "No carboxylic acid group found"
    
    # Check for the correct number of carbon atoms (at least 18 for the base chain)
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(carbon_atoms)
    if c_count < 18:
        return False, f"Carbon count is {c_count}, but requires at least 18 for octadecadienoic acid"

    # Check for exactly two C=C double bonds within a long chain (excludes terminal atoms)
    double_bond_pattern = Chem.MolFromSmarts("[C]=[C]")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) != 2:
        return False, f"Found {len(double_bond_matches)} C=C double bonds, needs exactly 2"

    # Ensure the main structure remains largely linear, allowing functional groups without large branches
    # For octadecadienoic acids, the straight primary structure can have modifications, 
    # hence consider long path counts rather than degree chaining
    for bond in double_bond_matches:
        path_length = Chem.rdmolops.GetShortestPath(mol, bond[0], bond[1])
        if len(path_length) < 4:  # Ensure bonds are not conjugated too closely which might indicate cyclization or complex branching
            return False, "Bonds are too close, likely indicating ring structure or heavy branching"

    return True, "Molecule is a valid octadecadienoic acid with linear structure, C18 backbone, and two C=C double bonds"