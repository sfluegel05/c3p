"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
from rdkit import Chem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    A 3alpha-hydroxy steroid is typically a steroid with a 3alpha-position hydroxyl group.

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
    
    # Define the steroid backbone pattern
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3(CCC4)C")
    
    # Check for steroid backbone presence
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Define the 3alpha-hydroxy pattern
    hydroxy_3alpha_pattern = Chem.MolFromSmarts("[C@@H]1(O)[C@@H](C)C2CCC3C4CCC(C3)C4(C)CC2C1")
    
    # Check for 3alpha-hydroxy group
    if not mol.HasSubstructMatch(hydroxy_3alpha_pattern):
        return False, "No 3alpha-hydroxy group found"
    
    return True, "Contains steroid backbone with 3alpha-position hydroxyl group"

# Examples of actual SMILES that should be classifiable
examples = [
    "[H][C@@]12C[C@H](O)CC[C@]1(C)[C@@]1([H])CC[C@]3(C)[C@]([H])(CC[C@@]3([H])[C@]1([H])[C@H](O)C2)[C@H](C)CCCC(C)CO", 
    "C1[C@@]2([C@]3(C[C@@H]([C@]4([C@]([C@@]3(CC[C@]2(C[C@@H](C1)O)[H])[H])(CC[C@@]4([C@@H](CCCC(C)(C)O)C)[H])[H])C)O)[H])C"
]

for ex in examples:
    print(f"SMILES: {ex} -> {is_3alpha_hydroxy_steroid(ex)}")