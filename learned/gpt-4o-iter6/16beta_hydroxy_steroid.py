"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    A 16beta-hydroxy steroid has a hydroxyl group at position 16 of the steroid backbone
    with a beta-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 16beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Generalized pattern for steroid backbone
    steroid_patterns = [
        Chem.MolFromSmarts("C1CCCC2C1CCC3C2CCC4C3(C)CC[C@@H]4C"),
        Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3CC[C@H]4C"),
        Chem.MolFromSmarts("C1CCCC2(C1)CCC3C2CCCC4C3CCC5=C4CCC=C5")
    ]
    
    # Ensure at least one steroid pattern matches
    if not any(mol.HasSubstructMatch(pattern) for pattern in steroid_patterns):
        return False, "No steroid backbone found"
    
    # Check for the 16beta-hydroxy group pattern
    beta_hydroxy_16_pattern = Chem.MolFromSmarts("[C@@H](O)[C@H]([CH2])[CH2][C@H]1CC[C@]2(C)[C@@H]([C@](O)([CH2])[CH2][C@@H]3[CH2][CH2][CH2][CH2]3)[C@@]1([CH3])[CH2][CH2]2")
    if not mol.HasSubstructMatch(beta_hydroxy_16_pattern):
        return False, "No 16beta-hydroxy group found"
    
    return True, "Contains a steroid backbone with 16beta-hydroxy group"

# Example usage
print(is_16beta_hydroxy_steroid("C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1C[C@H](O)[C@@H]2O"))