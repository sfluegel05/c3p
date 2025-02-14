"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: CHEBI:35692 dicarboxylic acid
"""

from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid is any carboxylic acid containing two carboxy groups (-COOH or -COO⁻),
    regardless of other functional groups or molecular weight.

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

    # Define carboxyl group pattern (both protonated and deprotonated forms)
    carboxyl_group_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H0-,OX2H1]")

    # Find all carboxyl groups in the molecule
    carboxyl_groups = mol.GetSubstructMatches(carboxyl_group_pattern)
    num_carboxyl_groups = len(carboxyl_groups)

    # Check if there are exactly two carboxyl groups
    if num_carboxyl_groups == 2:
        return True, "Molecule contains exactly two carboxyl groups"
    else:
        return False, f"Found {num_carboxyl_groups} carboxyl group(s), expected exactly 2"