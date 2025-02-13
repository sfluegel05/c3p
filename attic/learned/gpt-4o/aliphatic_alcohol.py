"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define pattern for aliphatic alcohols
    # CARBOXYLIC acids can also have carbonyls, check for non-aromatic carbon
    # Note: The pattern 'C[OH]' also matched diols (OH-C-CH(OH)-C) without aromaticity.
    pattern = Chem.MolFromSmarts("[$([CX4][OH]),$([CX3]=[CX3][OH]),$([CX2]#[CX0][OH]),$([C!r;R][OH])]")

    if pattern is None:
        return None, "Invalid SMARTS pattern"

    if mol.HasSubstructMatch(pattern):
        return True, "Molecule contains a hydroxyl group attached to an aliphatic carbon"

    return False, "No aliphatic alcohol pattern found"