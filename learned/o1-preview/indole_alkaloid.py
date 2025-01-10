"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: indole alkaloid
"""

from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid is an alkaloid containing an indole skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for indole ring
    indole_patterns = [
        'c1c2ccccc2[nH]c1',     # Indole with protonated nitrogen
        'c1c2ccccc2nc1',        # Indole with substituted nitrogen
        'c1nccc2ccccc12',       # Indole with nitrogen at different positions
        'c1nc2ccccc2c1',        # Another indole variation
        'n1c2ccccc2c[cH]c1',    # Indole with nitrogen at position 1
        'c1c2[nH]ccc2cc1',      # Indole with fused rings
        'c1c2nccc2cc1',         # Indole with nitrogen in different tautomer
        'c1c2cccnc2cc1',        # Indole with nitrogen in six-membered ring
    ]

    # Check if molecule contains any of the indole substructures
    for smarts in indole_patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            return True, "Contains indole skeleton"

    # If no patterns matched
    return False, "No indole skeleton found"

__metadata__ = {
    'chemical_class': {
        'id': '',  # CHEBI ID can be added if known
        'name': 'indole alkaloid',
        'definition': 'An alkaloid containing an indole skeleton.',
        'parents': []
    },
    # Additional metadata can be added here if needed
}