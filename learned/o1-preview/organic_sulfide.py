"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: CHEBI:15019 organic sulfide
"""
from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide based on its SMILES string.
    An organic sulfide is defined as a compound where a sulfur atom is bonded to two carbon atoms (R-S-R', R≠H).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an organic sulfide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the organic sulfide SMARTS pattern
    sulfide_pattern = Chem.MolFromSmarts("[#6]-[S;D2;H0]-[#6]")
    if not sulfide_pattern:
        return False, "Invalid SMARTS pattern"

    # Search for the sulfide pattern in the molecule
    if mol.HasSubstructMatch(sulfide_pattern):
        return True, "Contains an organic sulfide group (R-S-R', R≠H)"
    else:
        return False, "No organic sulfide group found"