"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol has a hydroxyl group (-OH) attached to a non-aromatic carbon,
    ideally in an open chain or a simple cycle, not involved in complex structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as an aliphatic alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for aliphatic alcohol
    # [C;R0] means non-aromatic, outside aromatic or complex ring systems
    aliphatic_oh_pattern = Chem.MolFromSmarts("[C&H2&H1][OX2H]")

    if mol.HasSubstructMatch(aliphatic_oh_pattern):
        return True, "Contains aliphatic hydroxyl group(s) fitting criteria for aliphatic alcohols"

    return False, "No fitting aliphatic hydroxyl group found for an aliphatic alcohol context"