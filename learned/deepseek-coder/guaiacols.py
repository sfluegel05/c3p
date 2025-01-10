"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: CHEBI:28417 guaiacols
Guaiacols are defined as any phenol carrying an additional methoxy substituent at the ortho-position.
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    A guaiacol is a phenol with a methoxy group at the ortho position.

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

    # Define a flexible SMARTS pattern for guaiacol: phenol with methoxy at ortho position
    # This pattern allows for additional substituents and modifications
    guaiacol_pattern = Chem.MolFromSmarts("[c;H0]([OH]):[c;H0]([OCH3])")
    
    # Check if the molecule matches the guaiacol pattern
    if mol.HasSubstructMatch(guaiacol_pattern):
        return True, "Contains a phenol with a methoxy group at the ortho position"
    else:
        return False, "Does not match the guaiacol pattern (phenol with methoxy at ortho position)"