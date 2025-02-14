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
    A primary alcohol is a compound in which a hydroxy group (-OH) is attached to a saturated carbon atom
    which has either three hydrogen atoms attached to it or only one other carbon atom and two hydrogen atoms attached to it.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a primary alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define primary alcohol SMARTS pattern
    primary_alcohol_pattern = Chem.MolFromSmarts("[#6;H2;X4][#8;H1]")  # Carbon with two hydrogens attached to OH group

    # Search for the pattern in the molecule
    if mol.HasSubstructMatch(primary_alcohol_pattern):
        return True, "Contains a primary alcohol group"
    else:
        return False, "No primary alcohol group found"