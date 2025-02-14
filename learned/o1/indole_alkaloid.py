"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: CHEBI:24842 indole alkaloid
"""

from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid is an alkaloid containing an indole skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for indole skeleton
    indole_smarts = 'c1cc2c([nH]c2cc1)'

    indole_pattern = Chem.MolFromSmarts(indole_smarts)

    if indole_pattern is None:
        return False, "Invalid indole SMARTS pattern"

    # Check if molecule contains indole skeleton
    if mol.HasSubstructMatch(indole_pattern):
        return True, "Contains indole skeleton (indole alkaloid)"
    else:
        return False, "No indole skeleton found"