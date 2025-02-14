"""
Classifies: CHEBI:47622 acetate ester
"""
from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule is an acetate ester based on its SMILES string.
    An acetate ester is any carboxylic ester where the carboxylic acid component is acetic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acetate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for an acetate ester group
    # This pattern matches a methyl group attached to a carbonyl carbon, which is single-bonded to an oxygen
    acetate_pattern = Chem.MolFromSmarts('[CH3][C](=O)O[*]')
    if acetate_pattern is None:
        return False, "Invalid SMARTS pattern for acetate ester"

    # Search for matches of the acetate ester pattern in the molecule
    matches = mol.GetSubstructMatches(acetate_pattern)
    if not matches:
        return False, "No acetate ester groups found"

    acetate_count = len(matches)

    return True, f"Contains {acetate_count} acetate ester group(s)"