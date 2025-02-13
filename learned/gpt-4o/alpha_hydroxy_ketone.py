"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
from rdkit import Chem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone contains a ketone group with a hydroxy group on the alpha-carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the ketone group (C=O)
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[#6]")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group found"

    # Look for hydroxyl group on the alpha carbon
    alpha_hydroxy_pattern = Chem.MolFromSmarts("[$([OH])[CX4;$([CX3](=O))]]")
    if not mol.HasSubstructMatch(alpha_hydroxy_pattern):
        return False, "No alpha-hydroxy group found relative to a ketone"

    return True, "Contains a ketone group with a hydroxy group on the alpha-carbon"