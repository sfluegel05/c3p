"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
from rdkit import Chem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    A 3beta-hydroxy steroid is characterized by a steroid structure with a hydroxyl
    group at the 3rd carbon in the beta orientation.

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
    
    # Define a simple pattern for the ABCD ring structure of steroids
    steroid_pattern = Chem.MolFromSmarts("C1CC2CCC3C(C(CCC4=CC(CCC34)C2)C1)C")
    if steroid_pattern is None:
        return None, "Failed to create a valid steroid pattern"

    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Define pattern for 3beta-hydroxy group
    # Match a tertiary alcohol pattern
    hydroxy_pattern = Chem.MolFromSmarts("C[C@@H](O)C")
    if hydroxy_pattern is None:
        return None, "Failed to create a valid hydroxy pattern"

    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3beta-hydroxy group found"

    return True, "Contains steroid backbone with 3beta-hydroxy group"