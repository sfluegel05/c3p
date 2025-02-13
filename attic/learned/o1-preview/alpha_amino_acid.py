"""
Classifies: CHEBI:33704 alpha-amino acid
"""
from rdkit import Chem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid is an amino acid in which the amino group is located on
    the carbon atom at the position alpha to the carboxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for alpha-amino acid
    # Pattern: amino group [N;!$(N=*);!$(N-[C]=O)] connected to alpha carbon [C;!$(C=*)],
    # which is connected to carboxylic acid [C(=O)[O;H]]
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[N;!$(N=*);!$(N-[C]=O)][C;!$(C=*)][C](=O)[O;H1]")

    if alpha_amino_acid_pattern is None:
        return False, "Invalid SMARTS pattern"

    match = mol.HasSubstructMatch(alpha_amino_acid_pattern)
    if match:
        return True, "Molecule matches alpha-amino acid pattern"
    else:
        return False, "Molecule does not match alpha-amino acid pattern"