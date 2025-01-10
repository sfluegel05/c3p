"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
from rdkit import Chem

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is defined as a polyunsaturated fatty acid with exactly three double bonds.

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

    # Check for carboxylic acid group: [CX3](=O)[OX1H]
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found; not a fatty acid"

    # Count carbon double bonds (C=C) ensuring they're part of a likely fatty acid structure
    double_bond_pattern = Chem.MolFromSmarts("C=CC=C")  # Ensure spacing denotes chain-like structure
    double_bond_matches = len(mol.GetSubstructMatches(double_bond_pattern))

    if double_bond_matches >= 3:
        return True, "Contains three or more suitably aligned double bonds; classified as trienoic fatty acid"
    else:
        return False, f"Has {double_bond_matches} properly aligned double bonds; not a trienoic fatty acid"