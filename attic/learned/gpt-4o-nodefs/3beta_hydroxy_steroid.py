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
    # This is a generalized pattern, since exact SMARTS is complex for this class; it matches steroid backbone and generic beta-hydroxylating
    steroid_pattern = Chem.MolFromSmarts("[C@@H]1([C@@H]2CC[C@]3(C)[C@@H](O)CC[C@]3(C)[C@@H]2[C@H](C)CC1)")

    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Does not match 3beta-hydroxy steroid pattern"

    return True, "Matches 3beta-hydroxy steroid pattern"