"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: Aliphatic Alcohol
"""
from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol is an alcohol derived from an aliphatic compound.
    This function checks if there is at least one hydroxyl group (-OH) attached to a non-aromatic, sp3-hybridized carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for hydroxyl group attached to non-aromatic, sp3-hybridized carbon
    alcohol_pattern = Chem.MolFromSmarts('[CX4;!Ar][OX2H]')

    # Find matches of the alcohol pattern
    matches = mol.GetSubstructMatches(alcohol_pattern)
    if matches:
        return True, "Contains aliphatic alcohol group(s)"
    else:
        return False, "No aliphatic alcohol groups found"