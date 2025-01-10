"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
from rdkit import Chem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    A 17beta-hydroxy steroid must have a hydroxy group at position 17 in beta-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Revised steroid pattern with more attention to backbone structure and flexibility
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CC3C4CCCCC4CC3C2')  # Typical tetracyclic structure of steroids
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Specialized pattern for the 17beta-hydroxy group
    hydroxy_17beta_pattern = Chem.MolFromSmarts('[C@H](O)[C@@H1C2CCC3C=CCC4C]=C(C)CC34C2C1]')  # Flexible matching for 17beta
    if not mol.HasSubstructMatch(hydroxy_17beta_pattern):
        return False, "No 17beta-hydroxy group found"

    return True, "Contains 17beta-hydroxy group with steroid backbone configuration"