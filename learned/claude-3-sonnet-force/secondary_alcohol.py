"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: CHEBI:35569 secondary alcohol

A secondary alcohol is a compound in which a hydroxy group, -OH, is attached to a saturated carbon atom which has two or more other carbon atoms attached to it.
"""

from rdkit import Chem
from rdkit.Chem import MolFromSmarts

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a secondary alcohol group, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for secondary alcohol group
    secondary_alcohol_pattern = MolFromSmarts("[CX4H2][OX2H]")

    # Check if the molecule has a substructure match for the secondary alcohol pattern
    if mol.HasSubstructMatch(secondary_alcohol_pattern):
        return True, "Contains a secondary alcohol group"

    return False, "No secondary alcohol group found"