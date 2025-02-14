"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: nitrile
"""
from rdkit import Chem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile is a compound having the structure R-C#N, where 'R' is a carbon substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a nitrile, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define nitrile SMARTS pattern:
    # [#6]-[#6]#[#7]: Carbon single-bonded to nitrile carbon, which is triple-bonded to nitrogen
    nitrile_pattern = Chem.MolFromSmarts("[#6]-[#6]#[#7]")
    if nitrile_pattern is None:
        return False, "Invalid nitrile SMARTS pattern"

    # Search for nitrile substructure
    if mol.HasSubstructMatch(nitrile_pattern):
        return True, "Contains nitrile group (R-C#N) per definition"
    else:
        return False, "Does not contain nitrile group (R-C#N) per definition"