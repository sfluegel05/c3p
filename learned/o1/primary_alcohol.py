"""
Classifies: CHEBI:15734 primary alcohol
"""
"""
Classifies: CHEBI:15734 primary alcohol
"""
from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.
    A primary alcohol is a compound in which a hydroxy group, -OH, is attached to a saturated carbon atom
    which has either three hydrogen atoms attached to it or only one other carbon atom and two hydrogen atoms attached to it.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for primary alcohol
    pattern = Chem.MolFromSmarts('[C;X4;H2,H3][O;H1]')

    # Search for substructure matches
    matches = mol.GetSubstructMatches(pattern)

    if len(matches) == 0:
        return False, "No primary alcohol group found"

    return True, f"Found {len(matches)} primary alcohol group(s)"