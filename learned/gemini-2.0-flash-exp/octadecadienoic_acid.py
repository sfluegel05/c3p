"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is a straight-chain C18 fatty acid with two C=C double bonds and a carboxylic acid group.

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
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for exactly two C=C double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) != 2:
         return False, f"Found {len(double_bond_matches)} double bonds, need exactly 2"

    # Check for 18 carbon atoms
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) != 18:
        return False, f"Found {len(carbon_atoms)} carbon atoms, need exactly 18"
    
    # Check that every carbon has at most 2 carbon neighbours
    for atom in carbon_atoms:
      carbon_neighbours = [neighbour for neighbour in atom.GetNeighbors() if neighbour.GetAtomicNum() == 6]
      if len(carbon_neighbours) > 2:
        return False, "Not a straight carbon chain"
    

    return True, "Contains 18 carbons in a straight chain, a carboxylic acid group and two double bonds"