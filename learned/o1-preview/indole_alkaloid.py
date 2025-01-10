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
    # Indole with protonated nitrogen
    indole_nh_smarts = 'c1c[nH]c2ccccc12'
    # Indole with substituted nitrogen
    indole_n_smarts = 'c1cnc2ccccc12'

    # Create RDKit molecule objects for SMARTS patterns
    indole_nh = Chem.MolFromSmarts(indole_nh_smarts)
    indole_n = Chem.MolFromSmarts(indole_n_smarts)

    # Check if molecule contains indole substructure
    if mol.HasSubstructMatch(indole_nh) or mol.HasSubstructMatch(indole_n):
        return True, "Contains indole skeleton"
    else:
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