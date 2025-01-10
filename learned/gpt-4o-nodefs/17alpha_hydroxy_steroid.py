"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
from rdkit import Chem

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason or structural match explanation
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved steroid core pattern
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CCC4C3(C)CCC4')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for 17-alpha hydroxyl group
    # Enhanced detection pattern
    hydroxy_17alpha_pattern = Chem.MolFromSmarts('O[C@H](C)CC[C@H]1CCC2C1CC[C@H]3[C@H]2CCC4=CCCC34')
    if not mol.HasSubstructMatch(hydroxy_17alpha_pattern):
        return False, "No 17alpha-hydroxy group found"
    
    return True, "Contains steroid backbone with 17alpha-hydroxy group"

# Example usage:
# Examples given in the task can be used to check for classification validity.
# result, reason = is_17alpha_hydroxy_steroid("your_example_smiles_here")
# print(result, reason)