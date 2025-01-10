"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
from rdkit import Chem

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    A 17alpha-hydroxy steroid has a hydroxyl group at the 17-alpha position on the steroid backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the steroid backbone SMARTS pattern
    # This pattern includes the four fused rings (phenanthrene skeleton) common to steroids
    steroid_backbone_smarts = "C1CCC2C(C1)CCC3C2CCC4(C3CCC4)"

    # Define SMARTS pattern for 17alpha-hydroxy group
    # The hydroxyl group should be attached to the D ring, which is the part of the molecule represented
    # by a cyclohexane attached to cyclopentane in steroids, at the 17th carbon.
    hydroxy_17alpha_smarts = "[C@]1([C@@H](O)CC2)CC3CC[C@H]4[C@@H](C2)CC[C@]34C1"

    # Check for the steroid backbone
    steroid_backbone = Chem.MolFromSmarts(steroid_backbone_smarts)
    if not mol.HasSubstructMatch(steroid_backbone):
        return False, "No steroid backbone found"

    # Check for presence of 17alpha-hydroxy group
    hydroxy_17alpha = Chem.MolFromSmarts(hydroxy_17alpha_smarts)
    if not mol.HasSubstructMatch(hydroxy_17alpha):
        return False, "No 17alpha-hydroxy group found"

    return True, "Contains steroid backbone with 17alpha-hydroxy group"