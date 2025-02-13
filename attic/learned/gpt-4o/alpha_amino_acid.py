"""
Classifies: CHEBI:33704 alpha-amino acid
"""
from rdkit import Chem

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

    # Improved SMARTS pattern for alpha-amino acids:
    # [N;H2] - Primary amine
    # [C@] - Ensures alpha carbon is chiral center (if specified)
    # [CX4] - Ensures alpha carbon is aliphatic (sp3)
    # [CX3](=O)[OX2H1] - Carboxylic acid group
    alpha_amino_acid_patterns = [
        "[N;H2]C[C;X4][C;X3](=O)[OX2H1]",                    # Basic structure
        "[N;H2][C@@H;!R][CX3](=O)[OX2H1]",                   # Chiral at alpha carbon
        "N[C@H](C(=O)O)C",                                   # Simple version with explicit structure
    ]

    for pattern in alpha_amino_acid_patterns:
        amino_acid_pattern = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(amino_acid_pattern):
            return True, "Molecule matches the structure of an alpha-amino acid"

    return False, "Molecule does not match the structure of an alpha-amino acid"