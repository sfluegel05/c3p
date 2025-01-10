"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
from rdkit import Chem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved steroid backbone pattern: 17 carbon atoms in 4 rings specific to steroids
    steroid_pattern = Chem.MolFromSmarts('C1CC[C@H]2C[C@H](O)CC[C@@]2([H])[C@H]1C3=CC4=CC=C[C@]3([H])[C@@]4([H])[C@H](O)2')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Basic steroid backbone not found"

    # Improved 11beta-hydroxy group pattern
    # This pattern now includes stereochemistry
    hydroxy_11beta_pattern = Chem.MolFromSmarts('[C@H]1(CC[C@H]2C[C@H](O)CC[C@@]2([H])[C@H]1C3)(O)[C@@]3(C)')
    if not mol.HasSubstructMatch(hydroxy_11beta_pattern):
        return False, "No 11beta-hydroxy group found with correct configuration"

    return True, "Contains basic steroid backbone with an 11beta-hydroxy group of the correct configuration"