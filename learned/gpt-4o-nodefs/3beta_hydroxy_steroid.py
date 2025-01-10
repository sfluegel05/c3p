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

    # Define SMARTS pattern for a steroid backbone with a 3beta-hydroxyl group
    # This is an attempt to capture the cyclopentanoperhydrophenanthrene structure with variability in stereochemistry and hydroxyl position at position 3.
    steroid_pattern = Chem.MolFromSmarts("[C3;R][C3;R][C3;R][C3;R][C;R]1[C;R][C;R](C)[C;R]([C;R]2[C;R]1[C;R](O)[C;R](=O)[C;R](C2)[C;R])")

    if mol.HasSubstructMatch(steroid_pattern):
        return True, "Matches 3beta-hydroxy steroid pattern"

    return False, "Does not match 3beta-hydroxy steroid pattern"