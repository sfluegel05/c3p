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

    # Look for steroid backbone using a generic steroid SMARTS pattern
    steroid_pattern = Chem.MolFromSmarts('C1C[C@H]2CC[C@@H]3C=CC(=O)CC[C@]3(C)C[C@@H]2[C@H]1C')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Look for 17beta-hydroxy group
    hydroxy_17beta_pattern = Chem.MolFromSmarts('[C@@H]1(O)[C@H](C)CCC[C@@H]2CCC[C@]12C')
    if not mol.HasSubstructMatch(hydroxy_17beta_pattern):
        return False, "No 17beta-hydroxy group pattern found"

    return True, "Contains 17beta-hydroxy group with correct steroid backbone configuration"