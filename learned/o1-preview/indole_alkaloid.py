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

    # Define indole SMARTS pattern (matches indole and N-substituted indoles)
    indole_smarts = 'c1nccc2ccccc12'  # Indole core structure, including N-substituted indoles
    indole = Chem.MolFromSmarts(indole_smarts)

    # Check if molecule has indole substructure
    if mol.HasSubstructMatch(indole):
        return True, "Contains indole skeleton"
    else:
        return False, "No indole skeleton found"

__metadata__ = {'chemical_class': { 'id': '',  # CHEBI ID can be added if known
                             'name': 'indole alkaloid',
                             'definition': 'An alkaloid containing an indole skeleton.',
                             'parents': []},
    # Additional metadata can be added here if needed
}