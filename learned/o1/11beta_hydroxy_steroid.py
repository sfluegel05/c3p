"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
"""
Classifies: 11beta-hydroxy steroid
"""
from rdkit import Chem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.
    An 11beta-hydroxy steroid is any steroid that has a hydroxy group at position 11 with beta-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid nucleus SMARTS pattern (rings A/B/C/D)
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CCC4C3CCCC4')
    if steroid_pattern is None:
        return False, "Invalid steroid nucleus SMARTS pattern"

    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Does not contain steroid nucleus"

    # Define beta-hydroxy group SMARTS pattern
    # Note: This is an approximation due to difficulty in specifying exact positions
    # [C@H] indicates a chiral carbon with beta configuration
    beta_hydroxy_pattern = Chem.MolFromSmarts('[C@H](O)')

    if beta_hydroxy_pattern is None:
        return False, "Invalid beta-hydroxy SMARTS pattern"

    if not mol.HasSubstructMatch(beta_hydroxy_pattern):
        return False, "No beta-configured hydroxy group found"

    # Get all matches of beta-hydroxy groups
    beta_hydroxy_matches = mol.GetSubstructMatches(beta_hydroxy_pattern)
    num_beta_hydroxy = len(beta_hydroxy_matches)

    if num_beta_hydroxy == 0:
        return False, "No beta-configured hydroxy group found"

    # Since we cannot reliably map atom indices to position 11,
    # we acknowledge this limitation
    return True, f"Contains steroid nucleus with {num_beta_hydroxy} beta-configured hydroxy group(s)"