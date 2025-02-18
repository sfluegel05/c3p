"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
from rdkit import Chem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    A 3alpha-hydroxy steroid is a 3-hydroxy steroid where the 3-hydroxy substituent is in the alpha-position.

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

    # Generalized pattern for steroid backbone: four connected cyclohexane rings.
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CC[C@@H]2C[C@H](C)CC[C@]2(C)[C@@H]1")  
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No steroid backbone found"
        
    # SMARTS pattern for 3alpha-hydroxy group in a chiral center
    # The @ indicates the stereochemistry (alpha)
    hydroxy_pattern = Chem.MolFromSmarts("[C@H]([O])C")  # Hydroxyl at a chiral center with alpha position
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3alpha-hydroxy group found"

    return True, "Contains steroid backbone with 3alpha-hydroxy group"

# Example usage
smiles_example = "C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(CC[C@]2(C[C@@H]([C@@H]1O)O)[H])[H])(CC[C@@]4([C@@H]([C@H]([C@@H](CC(C)C)O)O)C)[H])[H])C)[H])C"
print(is_3alpha_hydroxy_steroid(smiles_example))  # Expected: True with reason