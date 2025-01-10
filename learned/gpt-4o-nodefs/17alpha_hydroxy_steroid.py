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

    # Flexible steroid core pattern
    steroid_pattern = Chem.MolFromSmarts('[#6]1[#6][#6]2[#6](=O)[#6]([#6][#6]3[#6]([#6]1)[#6][#6][#6]4[#6]3(CC)CC[C@]4([#8])C)C2')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No flexible steroid backbone found"

    # Check for 17-alpha hydroxyl group
    hydroxy_17alpha_pattern = Chem.MolFromSmarts('O[C@H]1CCC2C1CC[C@H]3[C@H]2CCC4=CC(=O)CC[C@]34C') 
    if not mol.HasSubstructMatch(hydroxy_17alpha_pattern):
        return False, "No 17alpha-hydroxy group found"
    
    return True, "Contains flexible steroid backbone with 17alpha-hydroxy group"

# Example usage:
# Examples given in the task can be used to check for classification validity.
# result, reason = is_17alpha_hydroxy_steroid("your_example_smiles_here")
# print(result, reason)