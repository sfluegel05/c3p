"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
from rdkit import Chem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.

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

    # SMARTS pattern for detecting a general steroid backbone (tetracyclic core)
    steroid_pattern = Chem.MolFromSmarts("C1CC[C@H]2[C@H]1CC[C@H]3[C@@H]2CCC4=C3C[C@@H](O)CC4")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone detected"

    # SMARTS pattern for detecting 3alpha-hydroxy group specifically
    alpha_hydroxy_pattern = Chem.MolFromSmarts("[C@H](O)[C@@H]1CCC[C@H]2[C@H]1CCC3=C2CC[C@H](O)C3")
    if not mol.HasSubstructMatch(alpha_hydroxy_pattern):
        return False, "No 3alpha-hydroxy group detected"
    
    return True, "3alpha-hydroxy steroid structure identified"

# Example usage
smiles_example = "[H][C@@]12C[C@H](O)CC[C@]1(C)[C@@]1([H])C[C@H](O)[C@]3(C)[C@]([H])(CC[C@@]3([H])[C@]1([H])[C@H](O)C2)[C@H](C)CCCC(CO)CO"
result, reason = is_3alpha_hydroxy_steroid(smiles_example)
print(f"Is 3alpha-hydroxy steroid: {result}, Reason: {reason}")