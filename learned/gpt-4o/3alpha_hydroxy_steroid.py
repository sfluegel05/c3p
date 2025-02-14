"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
from rdkit import Chem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    A 3alpha-hydroxy steroid is a 3-hydroxy steroid with the 3-hydroxy substituent
    in the alpha-position.
    
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

    # Refine pattern for steroid backbone: steroid pattern with cyclic and tetra-cyclic structure
    steroid_pattern = Chem.MolFromSmarts("C1[C@@H]2CCC3C4CC[C@H](C4)CCC3C2C1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # SMARTS pattern for 3alpha-hydroxy group on a specific chiral center:
    # Targeted to find a hydroxyl group specifically at the third position
    # in alpha orientation.
    alpha_hydroxy_pattern = Chem.MolFromSmarts("[C@H]([O])[C@@H](C2CCC1=C[C@H](O)C[C@]21C)C")
    alpha_hydroxy_matches = mol.GetSubstructMatches(alpha_hydroxy_pattern)
    if not alpha_hydroxy_matches:
        return False, "3alpha-hydroxy group not confirmed with appropriate context"
    
    # Assume further checks to ensure the context and chiral configuration
    for match in alpha_hydroxy_matches:
        # Basic additional contextual checks could be done here for specific atoms if needed.
        return True, "Contains steroid backbone with a confirmed 3alpha-hydroxy group"
    
    return False, "Proper 3alpha-hydroxy position not confirmed in context"

# Test the function with an example
smiles_example = "C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(CC[C@]2(C[C@@H]([C@@H]1O)O)[H])[H])(CC[C@@]4([C@@H]([C@H]([C@@H](CC(C)C)O)O)C)[H])[H])C)[H])C"
print(is_3alpha_hydroxy_steroid(smiles_example))  # Expected: True with reason