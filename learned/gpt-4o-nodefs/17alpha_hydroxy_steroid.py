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

    # Typical steroid core (more flexible)
    steroid_pattern = Chem.MolFromSmarts('C1CCC2(C3CCC4CC(C3)CCC4C2C1)')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for 17-alpha hydroxyl group
    # Typical pattern around 17-alpha-hydroxy may be broader to allow flexibility
    hydroxy_17alpha_pattern = Chem.MolFromSmarts('[C@](O)([C;R1])')
    if not mol.HasSubstructMatch(hydroxy_17alpha_pattern):
        return False, "No 17alpha-hydroxy group found"
    
    return True, "Contains steroid backbone with 17alpha-hydroxy group"

# Example usage:
# Examples given in the task can be used to check for classification validity.
# result, reason = is_17alpha_hydroxy_steroid("your_example_smiles_here")
# print(result, reason)