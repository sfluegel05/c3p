"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
from rdkit import Chem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    A 3alpha-hydroxy steroid should have a steroid backbone and a hydroxyl group at the 3alpha position.

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
    
    # Generalized steroid backbone pattern allowing for various stereochemistry
    # Steroid backbones often look like a complex polycyclic structure
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3C(CC1)CCC3C2")  # Simplified steroid ring structure
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Pattern for a hydroxyl group at the 3alpha position
    # Allow flexibility in stereochemistry and correct relative position on ring
    hydroxy_3alpha_pattern = Chem.MolFromSmarts("[C@@H](O)C1CC[C@H]2C1C3CCC4C2C3CCC4")  # Simplified relative position
    if not mol.HasSubstructMatch(hydroxy_3alpha_pattern):
        return False, "No 3alpha-hydroxy group found"
    
    return True, "Contains steroid backbone with 3alpha-position hydroxyl group"

# Testing with examples
examples = [
    "[H][C@@]12C[C@H](O)CC[C@]1(C)[C@@]1([H])CC[C@]3(C)[C@]([H])(CC[C@@]3([H])[C@]1([H])[C@H](O)C2)[C@H](C)CCCC(C)CO", 
    "C1[C@@]2([C@]3(C[C@@H]([C@]4([C@]([C@@]3(CC[C@]2(C[C@@H](C1)O)[H])[H])(CC[C@@]4([C@@H](CCCC(C)(C)O)C)[H])[H])C)O)[H])C"
]

for ex in examples:
    result, reason = is_3alpha_hydroxy_steroid(ex)
    print(f"SMILES: {ex} -> Result: {result}, Reason: {reason}")