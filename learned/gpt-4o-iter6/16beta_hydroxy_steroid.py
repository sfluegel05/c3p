"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
from rdkit import Chem

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

    # Correct pattern for a general steroid backbone using 6-6-6-5 ring system
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4C(C3)CCC4C2C1")

    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Pattern to match the hydroxyl group at position 16 with beta configuration.
    # Assuming that the pattern attaches correctly considering stereochemistry,
    beta_hydroxy_16_pattern = Chem.MolFromSmarts("C[C@@H](O)C(C)C1CC[C@]2(C)CC[C@H](C2)O")

    if not mol.HasSubstructMatch(beta_hydroxy_16_pattern):
        return False, "No 16beta-hydroxy group found"
    
    return True, "Contains a steroid backbone with 16beta-hydroxy group"

# Example usage
print(is_16beta_hydroxy_steroid("C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1C[C@H](O)[C@@H]2O"))