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
        bool: True if molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for ortho-methoxy group adjacent to phenol
    # This pattern searches for a methoxy and hydroxy group on the same benzene ring
    # where the methoxy is ortho with respect to the hydroxy
    guaiacol_pattern = Chem.MolFromSmarts("c1cc(OC)cc(O)c1")
    
    # Check for match of guaiacol pattern
    if mol.HasSubstructMatch(guaiacol_pattern):
        return True, "Contains ortho-methoxy group adjacent to phenol on an aromatic ring"
    
    return False, "The structure does not match the guaiacol pattern"