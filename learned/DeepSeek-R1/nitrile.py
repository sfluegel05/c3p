"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: CHEBI:16330 nitrile
"""
from rdkit import Chem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile is a compound with the structure RC#N where R is any organic substituent.

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
    
    # Define SMARTS pattern for nitrile group (RC#N where R is carbon-based)
    # The pattern matches a carbon connected to another carbon via single bond, then triple bond to N
    nitrile_pattern = Chem.MolFromSmarts("[C][C]#[N]")
    
    # Check for presence of nitrile group with carbon substituent
    if mol.HasSubstructMatch(nitrile_pattern):
        return True, "Contains a nitrile group (RC#N)"
    else:
        return False, "No nitrile group (RC#N) found"