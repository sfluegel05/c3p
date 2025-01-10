"""
Classifies: CHEBI:134251 guaiacols
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    A guaiacol is any phenol carrying an additional methoxy substituent at the ortho-position.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for guaiacols (phenol with ortho-methoxy group)
    guaiacol_pattern = Chem.MolFromSmarts('[OH]c1c([OCH3])cccc1')
    if guaiacol_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check if molecule matches the guaiacol pattern
    if mol.HasSubstructMatch(guaiacol_pattern):
        return True, "Contains phenol with ortho-methoxy group at ortho-position"
    else:
        return False, "Does not contain phenol with ortho-methoxy group at ortho-position"