"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: secondary alcohol
"""
from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.
    A secondary alcohol has a hydroxy group (-OH) attached to a saturated carbon atom
    which is connected to exactly two other carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for secondary alcohol
    # [C:1] - saturated carbon atom
    # [O:2] - hydroxyl oxygen
    # Attached to two other carbon atoms
    pattern = Chem.MolFromSmarts("[CX4;H1]([OX2H])([CX4H])([CX4H])")
    if pattern is None:
        return False, "Invalid SMARTS pattern"

    # Search for matches
    matches = mol.GetSubstructMatches(pattern)
    if matches:
        return True, "Contains a secondary alcohol functional group"
    else:
        return False, "No secondary alcohol functional group found"