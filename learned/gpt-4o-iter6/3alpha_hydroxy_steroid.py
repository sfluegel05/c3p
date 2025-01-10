"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
from rdkit import Chem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    A 3alpha-hydroxy steroid is defined as a steroid with a hydroxyl group at the 3-position
    in an alpha orientation.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Flexible steroid backbone pattern, allowing for various stereocenters and saturation.
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C(C)CC2)C1")  # Simplified representation of the steroid nucleus
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone detected"
    
    # Pattern specifically for 3alpha-hydroxy group.
    # The hydroxyl needs to be in the alpha position at the 3-carbon.
    three_alpha_hydroxy_pattern = Chem.MolFromSmarts("[C@@H](O)C(C)C")  # Indicates alpha OH on C3
    if not mol.HasSubstructMatch(three_alpha_hydroxy_pattern):
        return False, "No 3alpha-hydroxy group detected"
    
    return True, "3alpha-hydroxy steroid structure identified"

# Example usage
smiles_example = "[H][C@@]12C[C@H](O)CC[C@]1(C)[C@@]1([H])C[C@H](O)[C@]3(C)[C@]([H])(CC[C@@]3([H])[C@]1([H])[C@H](O)C2)[C@H](C)CCCC(CO)CO"
result, reason = is_3alpha_hydroxy_steroid(smiles_example)
print(f"Is 3alpha-hydroxy steroid: {result}, Reason: {reason}")