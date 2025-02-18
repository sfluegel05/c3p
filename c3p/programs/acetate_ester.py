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

    # Define the SMARTS pattern for an acetate ester
    # This pattern looks for a methyl group [CH3], attached to a carbonyl carbon [C](=O),
    # which is connected via an oxygen to any non-hydrogen atom O[!H0]
    acetate_pattern = Chem.MolFromSmarts("[CH3][C](=O)O[!H0]")
    if acetate_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Search for matches of the acetate ester pattern in the molecule
    matches = mol.GetSubstructMatches(acetate_pattern)
    if matches:
        num_acetates = len(matches)
        return True, f"Contains {num_acetates} acetate ester group(s)"
    else:
        return False, "No acetate ester groups found"