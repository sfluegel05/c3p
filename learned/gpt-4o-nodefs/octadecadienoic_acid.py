"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
from rdkit import Chem

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is a fatty acid with a chain of 18 carbon atoms and exactly 2 double bonds,
    ending with a carboxylic acid group.

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

    # Check for carboxylic acid group presence
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count the number of double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_count = len(mol.GetSubstructMatches(double_bond_pattern))
    if double_bond_count != 2:
        return False, f"Expected exactly 2 double bonds, found {double_bond_count}"

    # Count total number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 18:
        return False, f"Expected 18 carbon atoms, found {c_count}"
    
    return True, "Molecule is an octadecadienoic acid with 18 carbons and 2 double bonds"