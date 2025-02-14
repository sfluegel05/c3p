"""
Classifies: CHEBI:134251 guaiacols
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    Guaiacols are phenols with an additional methoxy substituent at the ortho-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for ortho-methoxyphenol (guaiacol structure)
    guaiacol_pattern = Chem.MolFromSmarts("c1cc(OC)c(O)c1")
    
    if mol.HasSubstructMatch(guaiacol_pattern):
        return True, "Contains ortho-methoxy group adjacent to phenol on an aromatic ring"

    return False, "The structure does not match the guaiacol pattern"