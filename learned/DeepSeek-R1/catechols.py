"""
Classifies: CHEBI:33566 catechols
"""
"""
Classifies: CHEBI:33853 catechol
"""
from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule is a catechol (contains an o-diphenol group) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains adjacent aromatic hydroxyl groups, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for o-diphenol (two adjacent aromatic hydroxyl groups)
    # Pattern matches two hydroxyl groups attached to adjacent aromatic carbons
    catechol_pattern = Chem.MolFromSmarts("[OH]c:c[OH]")
    
    # Check for the presence of the pattern
    if mol.HasSubstructMatch(catechol_pattern):
        return True, "Contains adjacent aromatic hydroxyl groups (o-diphenol)"
    else:
        return False, "No o-diphenol component found"