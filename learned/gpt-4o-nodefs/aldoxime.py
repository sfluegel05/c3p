"""
Classifies: CHEBI:22307 aldoxime
"""
from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime contains a carbon-nitrogen double bond (C=N) and a hydroxyl group
    directly attached to the nitrogen (C=N-O structure), where the carbon is near-terminal
    mimicking an aldehyde structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aldoxime, False otherwise
        str: Reason for classification
    """

    # Parse SMILES into RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more specific SMARTS pattern for an aldoxime
    # The carbon should be terminal or near-terminal, akin to how it's present in an aldehyde.
    aldoxime_pattern = Chem.MolFromSmarts("[CX3H1]=[NX2]O")  # Restricting ensuring a single hydrogen bonded to enable it resembles aldehyde C.

    if mol.HasSubstructMatch(aldoxime_pattern):
        return True, "Contains the characteristic C=N-O structure of aldoximes with aldehyde-like carbon"

    return False, "Does not contain the characteristic C=N-O structure of aldoximes"