"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester is the ester derivative of an alpha-amino acid, formed by the condensation with an alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern: Alpha-amino acid with ester linkage
    alpha_amino_pattern = "[CX4H2,CX3H](-[NH2,NH])[CX3](=O)[O][C]"  # Alpha carbon, amino group, and ester linkage

    alpha_amino_mol = Chem.MolFromSmarts(alpha_amino_pattern)

    if alpha_amino_mol is None:
        return None, None  # Error in pattern creation

    if mol.HasSubstructMatch(alpha_amino_mol):
        return True, "Contains an alpha-amino acid backbone with ester linkage"

    return False, "Does not match the alpha-amino acid ester pattern"