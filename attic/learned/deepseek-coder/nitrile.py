"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: CHEBI:18379 nitrile
"""
from rdkit import Chem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile is a compound with the structure RC#N, where R is any organic group.

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

    # Look for the nitrile functional group (C#N)
    nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
    if not mol.HasSubstructMatch(nitrile_pattern):
        return False, "No nitrile functional group (C#N) found"

    return True, "Contains a nitrile functional group (C#N)"