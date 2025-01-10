"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
from rdkit import Chem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more appropriate SMARTS pattern for a steroid backbone 
    # with a flexible but definitive hydroxyl group on the third carbon in a beta orientation
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3CC[C@H](C4=CC(=O)CC[C@]24C)C3C1O")

    if mol.HasSubstructMatch(steroid_pattern):
        return True, "Matches 3beta-hydroxy steroid pattern"

    return False, "Does not match 3beta-hydroxy steroid pattern"