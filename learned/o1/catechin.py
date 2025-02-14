"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    A catechin is a type of flavonoid with a flavan-3-ol skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a catechin, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the flavan-3-ol skeleton without stereochemistry
    flavan3ol_smarts = """
    [$([#6]1:c:[cH]:c:[cH]:c1)-!@[#6]2-[#6]-[#8]-c3c([#6]2)-c:c:c:c3]-[#8]
    """
    flavan3ol_pattern = Chem.MolFromSmarts(flavan3ol_smarts.strip())
    if flavan3ol_pattern is None:
        return False, "Error in SMARTS pattern"

    # Check if the molecule contains the flavan-3-ol skeleton
    if mol.HasSubstructMatch(flavan3ol_pattern):
        return True, "Contains flavan-3-ol skeleton characteristic of catechins"
    else:
        return False, "Flavan-3-ol skeleton not found"