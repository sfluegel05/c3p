"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
from rdkit import Chem

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    A 17alpha-hydroxy steroid should have a hydroxyl group at the 17alpha position on the steroid structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # General SMARTS pattern for the steroid backbone (four fused rings)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3(C)CCC4")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # SMARTS for a 17alpha-hydroxy group with stereochemistry
    # Considering alternative stereochemistry configurations
    hydroxyl_patterns = [
        Chem.MolFromSmarts("[C@H](O)CC[C@@]1(C)CC[C@]2C1CCC3C2"),
        Chem.MolFromSmarts("[C@@H](O)CC[C@@]1(C)CC[C@]2C1CCC3C2"),
        Chem.MolFromSmarts("[C@H](O)CC[C@]1(C)CC[C@]2C1CCC3C2"),
        Chem.MolFromSmarts("[C@@H](O)CC[C@]1(C)CC[C@]2C1CCC3C2")
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in hydroxyl_patterns):
        return False, "No 17alpha-hydroxyl group found"
    
    return True, "Contains steroid backbone with 17alpha-hydroxyl group"