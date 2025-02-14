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
    A nitrile contains a carbon triple-bonded to a nitrogen atom (C#N).

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

    # Define nitrile SMARTS pattern (carbon triple-bonded to nitrogen)
    nitrile_pattern = Chem.MolFromSmarts("[CX2]#N")
    if nitrile_pattern is None:
        return False, "Invalid nitrile SMARTS pattern"
    
    # Search for nitrile substructure
    if mol.HasSubstructMatch(nitrile_pattern):
        return True, "Contains nitrile functional group (C#N)"
    else:
        return False, "Does not contain nitrile functional group (C#N)"