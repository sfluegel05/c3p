"""
Classifies: CHEBI:33566 catechols
"""
from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule is a catechol based on its SMILES string.
    A catechol contains an o-diphenol component (two hydroxyl groups on adjacent carbons on an aromatic ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a catechol moiety, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More precise SMARTS pattern for ortho-diphenol
    # Focusing on benzene-like systems due to common occurrence in catechols
    catechol_pattern = Chem.MolFromSmarts("c1cc(O)c(O)c1")  # Specifically matches benzene rings with ortho-hydroxy groups

    # Check the pattern
    if mol.HasSubstructMatch(catechol_pattern):
        return True, "Contains o-diphenol component in an aromatic ring system"
    
    return False, "Does not contain o-diphenol component in an aromatic ring system"