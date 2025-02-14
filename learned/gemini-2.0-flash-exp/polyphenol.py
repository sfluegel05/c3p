"""
Classifies: CHEBI:26195 polyphenol
"""
from rdkit import Chem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.
    A polyphenol has 2 or more benzene rings, each with at least one directly attached -OH group.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a polyphenol, False otherwise.
        str: Reason for classification.
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a benzene ring with a directly attached hydroxyl group
    polyphenol_pattern = Chem.MolFromSmarts("c1([c](O)ccccc1)")

    # Count the number of times the pattern is found
    matches = mol.GetSubstructMatches(polyphenol_pattern)
    num_matches = len(matches)

    # Classify based on the number of matches
    if num_matches >= 2:
        return True, "Contains two or more benzene rings, each with at least one directly attached -OH group."
    else:
      return False, f"Found {num_matches} benzene rings with directly attached -OH groups, needs at least two."