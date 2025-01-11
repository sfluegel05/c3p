"""
Classifies: CHEBI:18379 nitrile
"""
from rdkit import Chem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile is characterized by the presence of a C#N group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrile, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define nitrile group pattern (C triple bonded to N)
    nitrile_pattern = Chem.MolFromSmarts("C#N")
    if mol.HasSubstructMatch(nitrile_pattern):
        return True, "Contains nitrile group (C#N)"
    else:
        return False, "No nitrile group found"