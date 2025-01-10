"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester is defined as having an alpha-carbon connected to an amino group
    and a carboxyl group esterified with an alcohol.

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

    # Pattern: Alpha carbon connected to amino and esterified carboxyl group
    # SMARTS: [NX3;H2,H1,H0][CX4][C](=O)O[!H] (captures ester linkage with an alcohol)
    amino_group = "[NX3;H2,H1,H0]"  # Amino group, allowing for substitutions
    alpha_carbon = "[CX4]"           # Tetrahedral alpha-carbon
    ester_linkage = "[C](=O)O[!H]"   # Ester linkage indicating esterified carboxyl group
    alpha_amino_acid_ester_pattern = Chem.MolFromSmarts(amino_group + alpha_carbon + ester_linkage)
    
    if not alpha_amino_acid_ester_pattern:
        return None, None  # Error in pattern creation

    if mol.HasSubstructMatch(alpha_amino_acid_ester_pattern):
        return True, "Contains an alpha-amino acid backbone with ester linkage"
    
    return False, "Does not match the alpha-amino acid ester pattern"