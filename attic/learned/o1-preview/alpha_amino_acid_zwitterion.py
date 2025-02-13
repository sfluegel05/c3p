"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    An alpha-amino-acid zwitterion is characterized by a central alpha carbon connected to a protonated amino group ([NH3+]) and a deprotonated carboxylate group ([C(=O)[O-]]).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino-acid zwitterion, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for an alpha-amino-acid zwitterion
    pattern = Chem.MolFromSmarts('[N+;H3][C;H]([*!H])[C](=O)[O-]')
    if pattern is None:
        return False, "Invalid SMARTS pattern"

    if mol.HasSubstructMatch(pattern):
        return True, "Contains alpha-amino-acid zwitterion motif"
    else:
        return False, "Missing alpha-amino-acid zwitterion motif"