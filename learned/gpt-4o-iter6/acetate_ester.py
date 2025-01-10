"""
Classifies: CHEBI:47622 acetate ester
"""
from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.
    An acetate ester specifically contains the esterified acetic acid moiety (CH3COO-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an acetate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a refined acetate ester pattern: methyl group as part of the ester
    acetate_pattern = Chem.MolFromSmarts("[CH3]C(=O)O[!$([#6](=O)O)]")

    # Check if molecule contains the refined acetate ester pattern
    if mol.HasSubstructMatch(acetate_pattern):
        return True, "Contains acetic acid ester group (CH3COO-)"

    return False, "Does not contain the acetic acid ester group (CH3COO-)"

# Testing examples will proceed similarly with multiple SMILES inputs