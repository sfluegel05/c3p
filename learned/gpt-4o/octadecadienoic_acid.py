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

    # Check for carboxylic acid group at one end
    carboxylic_acid_group = Chem.MolFromSmarts("C(=O)O")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid_group)
    if len(carboxylic_matches) < 1:
        return False, "No carboxylic acid group found"
    
    # Check for exactly 18 carbon atoms in total
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) != 18:
        return False, f"Carbon count is {len(carbon_atoms)}, but requires exactly 18 for octadecadienoic acid"

    # Check for exactly two C=C double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) != 2:
        return False, f"Found {len(double_bond_matches)} C=C double bonds, needs exactly 2"

    # Ensure that the main carbon chain remains linear
    # Check the connectivity where each C has only one or two heavy atom neighbors
    chain_pattern = Chem.MolFromSmarts("[C;!R](=O)([!#1])[CH2]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Non-linear primary carbon chain or excessive branching found"

    # Check double-bond separation - ensure there's a minimum number of carbons (3) separating the two C=C 
    # to prevent conjugation misinterpretation or close packing
    for bond_a, bond_b in double_bond_matches:
        path_length = len(Chem.rdmolops.GetShortestPath(mol, bond_a, bond_b))
        if path_length <= 3:
            return False, "Double bonds are too close, could indicate incorrect structure interpretation"

    return True, "Molecule is a valid octadecadienoic acid with C18 backbone and two C=C double bonds"