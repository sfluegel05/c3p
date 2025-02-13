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
    
    # Define SMARTS pattern for the steroid backbone (tetracyclic structure)
    steroid_backbone_smarts = "C1CCC2C3CCC4[C@H]3CC[C@@]4([H])C2[C@@H]1"
    
    # Define SMARTS pattern for the 17alpha-hydroxy group completely
    hydroxy_17alpha_smarts = "C1CCC2C3CCC4[C@H]3CC[C@@]4(O)[C@H](C(=O)O)C2[C@@H]1"

    # Check for steroid backbone
    steroid_backbone = Chem.MolFromSmarts(steroid_backbone_smarts)
    if not mol.HasSubstructMatch(steroid_backbone):
        return False, "No steroid backbone found"
    
    # Check for presence of 17alpha-hydroxy group
    hydroxy_17alpha = Chem.MolFromSmarts(hydroxy_17alpha_smarts)
    if not mol.HasSubstructMatch(hydroxy_17alpha):
        return False, "No 17alpha-hydroxy group found"

    return True, "Contains steroid backbone with 17alpha-hydroxy group"