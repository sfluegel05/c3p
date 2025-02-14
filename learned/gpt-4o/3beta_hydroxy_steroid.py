"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
from rdkit import Chem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    A 3beta-hydroxy steroid should have the stereo-specific hydroxyl group at carbon 3 in the beta position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3beta-hydroxy steroid, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible steroid backbone pattern
    steroid_pattern = Chem.MolFromSmarts('C1CC[C@H]2C[C@H]3[C@@H](CC2)C[C@H]4CC[C@@H](C1)C34')
    
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
        
    # Define the 3beta-hydroxy pattern with beta orientation at C3 explicitly
    hydroxy_3beta_pattern = Chem.MolFromSmarts('[C@H]1C(C2CC[C@H]3[C@@H](C2)CC[C@@H]4CC[C@H](C1)C34)[OH]')
    if not mol.HasSubstructMatch(hydroxy_3beta_pattern):
        return False, "3beta-hydroxy group not found"
    
    return True, "Contains steroid backbone and a 3beta-hydroxy group"