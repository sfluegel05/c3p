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

    # Define a general steroid backbone SMARTS pattern
    # This pattern captures the four rings with stereochemistry flexibility
    steroid_backbone_smarts = "[#6]1[#6](=[#8])-[#6]2-[#6](=[#8])-[#6]3-[#6]4[#6](=[#8])-[#6]([#6]2)[#6]([#6]3)[#6](=[#8])[#6](1)[#6](=[#8])[#6]4"
    
    # Define SMARTS pattern for a 17alpha-hydroxy group
    # This aims to check for an OH group on the D ring with specific stereochemistry
    hydroxy_17alpha_smarts = "[C@@]1([C@@H](O)[C@H]2)[C@H]3CC[C@@H]4[C@@H]([C@H]2)CC[C@]34C1"
    
    # Check for the steroid backbone
    steroid_backbone = Chem.MolFromSmarts(steroid_backbone_smarts)
    if not mol.HasSubstructMatch(steroid_backbone):
        return False, "No steroid backbone found"

    # Check for presence of 17alpha-hydroxy group
    hydroxy_17alpha = Chem.MolFromSmarts(hydroxy_17alpha_smarts)
    if not mol.HasSubstructMatch(hydroxy_17alpha):
        return False, "No 17alpha-hydroxy group found"

    return True, "Contains steroid backbone with 17alpha-hydroxy group"