"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid is an amino acid where the amino group is acetylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an N-acetyl-amino acid, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for N-acetyl-amino acid
    # Pattern breakdown:
    # [CH3] - methyl group of acetyl
    # [C](=O) - carbonyl carbon of acetyl group
    # [N] - nitrogen atom
    # [C] - alpha carbon
    # [C](=O)[O,H1,H0] - carboxyl group (handles carboxylic acid or carboxylate forms)
    pattern = Chem.MolFromSmarts('[CH3][C](=O)[N][C][C](=O)[O,H1,H0]')

    if mol.HasSubstructMatch(pattern):
        return True, "Molecule matches N-acetyl-amino acid pattern"
    else:
        return False, "Molecule does not match N-acetyl-amino acid pattern"