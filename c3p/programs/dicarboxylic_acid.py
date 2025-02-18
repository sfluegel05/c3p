"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: CHEBI:35692 dicarboxylic acid
"""

from rdkit import Chem
from rdkit.Chem import rdqueries

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid is any carboxylic acid containing two carboxy groups (-COOH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dicarboxylic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define carboxyl group pattern (-C(=O)O[H or -])
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-1]")
    if carboxyl_pattern is None:
        return False, "Failed to create carboxyl group pattern"

    # Find all carboxyl groups in the molecule
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    num_carboxyl_groups = len(carboxyl_matches)

    # Check if there are exactly two carboxyl groups
    if num_carboxyl_groups == 2:
        return True, "Molecule contains exactly two carboxyl groups"
    else:
        return False, f"Found {num_carboxyl_groups} carboxyl group(s), expected exactly 2"