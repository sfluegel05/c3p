"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: 3β-hydroxy steroid
Definition: A 3-hydroxy steroid in which the 3-hydroxy substituent is in the β‐position.
"""

from rdkit import Chem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3β-hydroxy steroid based on its SMILES string.
    The classifier verifies:
      1. The presence of a steroid nucleus (i.e. a fused tetracyclic ring system).
      2. The presence of at least one β-oriented hydroxy group, indicated by the [C@@H](O) pattern.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 3β-hydroxy steroid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string to create a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the steroid nucleus.
    # This pattern represents a common steroid (gonane) core,
    # although it may not catch all steroid variants.
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C1CC[C@H]3[C@@H]2CC[C@]4(C)[C@@H]3CC[C@]4(C)")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Steroid nucleus not found"
    
    # Define a SMARTS pattern for a beta-oriented hydroxyl group.
    # In many SMILES, a beta OH is encoded as attached to a stereocenter using [C@@H](O).
    beta_hydroxy_pattern = Chem.MolFromSmarts("[C@@H](O)")
    if not mol.HasSubstructMatch(beta_hydroxy_pattern):
        return False, "Beta-oriented hydroxyl group not found"
    
    return True, "Molecule contains a steroid nucleus and a beta-oriented (3β) hydroxyl group"

# Example usage (uncomment to test):
# test_smiles = "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CC[C@@H](O)C(C)C"  # (24R)-24-hydroxycholesterol
# result, reason = is_3beta_hydroxy_steroid(test_smiles)
# print(result, reason)