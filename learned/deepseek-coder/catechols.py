"""
Classifies: CHEBI:33566 catechols
"""
"""
Classifies: CHEBI:33567 catechol
"""
from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule is a catechol based on its SMILES string.
    A catechol is any compound containing an o-diphenol component (two hydroxyl groups on adjacent carbons of a benzene ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the catechol substructure pattern (o-diphenol)
    catechol_pattern = Chem.MolFromSmarts("[c;H1]1:[c;H1]([OH]):[c;H1]:[c;H1]:[c;H1]:[c;H1]:1[OH]")
    
    # Check if the molecule contains the catechol pattern
    if mol.HasSubstructMatch(catechol_pattern):
        return True, "Contains an o-diphenol component (catechol)"
    else:
        return False, "No o-diphenol component found"