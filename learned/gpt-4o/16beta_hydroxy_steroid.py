"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
from rdkit import Chem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    A 16beta-hydroxy steroid is defined as a steroid with an -OH group at position 16 in beta-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 16beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # General steroid backbone pattern: four fused rings
    steroid_backbone_pattern = Chem.MolFromSmarts("[C&R1]1[C&R2][C&R2][C&R3][C&R3][C&R4][C&R4][C&R5]2[C&R5][C&R6][C&6]([C&R6][C&R7][C&R7][C&R8][C&R8]1)CCC2")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No general steroid backbone found"
    
    # Beta-OH group at position 16 using the correct stereochemistry
    # Assuming that the 16th carbon is an alpha or beta position:
    beta_oh_16_pattern = Chem.MolFromSmarts("[C@H](O)[C@]([C&R1])([C&R2][C&R3][C&R3][C&R4][C&R4]C)")
    if mol.HasSubstructMatch(beta_oh_16_pattern):
        return True, "Contains 16beta-hydroxy group with correct steroid backbone"
    
    return False, "Does not match 16beta-hydroxy steroid criteria"