"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
from rdkit import Chem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string,
    requiring a steroid backbone and a stereospecific 11beta-hydroxy group.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more comprehensive steroid backbone pattern
    # A basic A/B/C/D ring fusion structure (assume [6-6-6-5] with flexibility)
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CC2CCC3C[C@@H](O)CCC3C2C1")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "Steroid backbone not found"

    # Check for a specific 11beta-hydroxy configuration (assuming it to be more flexible)
    # Use SMARTS pattern that specifies the exact configuration
    hydroxy_11beta_pattern = Chem.MolFromSmarts("[C@H](O)[C@@H]1CCC2=C1CCC3C2CCC4C3CCC(C4)C")
    if not mol.HasSubstructMatch(hydroxy_11beta_pattern):
        return False, "No 11beta-hydroxy group with correct configuration found"

    return True, "Contains steroid backbone with an 11beta-hydroxy group of the expected configuration"