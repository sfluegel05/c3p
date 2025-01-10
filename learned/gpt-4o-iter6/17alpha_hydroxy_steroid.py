"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
from rdkit import Chem

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    A 17alpha-hydroxy steroid has a hydroxyl group at the 17alpha position on the steroid structure.

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
    
    # Improved SMARTS pattern for steroid backbone, accounting for stereochemistry
    # Note: "[#6]" denotes a carbon atom to cover cases where stereochemistry may cause issues
    steroid_pattern = Chem.MolFromSmarts('C1CC[C@]2(C)C(C)([C@@H]([C@@H]2CC1)[C@H](O))C')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Improved SMARTS pattern for the 17alpha-hydroxy group
    # Often involves "[C@]1(O)[C@H](C)"
    hydroxyl_pattern = Chem.MolFromSmarts('[C@H](O)[C@H]1')
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No 17alpha-hydroxyl group found"
    
    # Additional checks or patterns can be added to ensure stereochemistry is captured correctly
    
    return True, "Contains steroid backbone with 17alpha-hydroxyl group"

# Example usage with a known 17alpha-hydroxy steroid
smiles_example = '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(=O)CO'
print(is_17alpha_hydroxy_steroid(smiles_example))