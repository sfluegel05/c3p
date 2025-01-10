"""
Classifies: CHEBI:33704 alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid has an amino group on the carbon atom at the position alpha to the carboxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for main features of an alpha amino acid
    # [NX3;H2,H1;!$(NC=O)] - Looks for an amine group not already part of an amide linkage
    # [CX4H1,O] - Ensures alpha carbon is sp3 hybridized and linked appropriately (can be complex pattern)
    # [CX3](=O)[OX2H1] - Looks for a carboxylic acid group
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]C[C;X4][CX3](=O)[OX2H1]")

    # Check if the molecule matches the alpha amino acid pattern
    if mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return True, "Molecule has the structure of an alpha-amino acid"

    return False, "Molecule does not match the pattern of an alpha-amino acid"