"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    A 16beta-hydroxy steroid is a steroid with an -OH group at position 16 in beta-configuration.

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
    
    # Steroid backbone pattern: look for four fused rings typical in steroids
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CC3C4CCC(C3)C(C4)C2")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No general steroid backbone found"
    
    # Look for beta-OH at position 16
    beta_oh_16_pattern = Chem.MolFromSmarts("[C@H](O)[C@@H](C)[CH2][CH2][CH2][CH2]C")
    if not mol.HasSubstructMatch(beta_oh_16_pattern):
        return False, "No beta-OH at position 16 found"
    
    # Verify the OH group at the beta position
    if mol.HasSubstructMatch(beta_oh_16_pattern):
        return True, "Contains beta-hydroxy group on steroid skeleton at position 16 in beta-configuration"

    return False, "Does not match 16beta-hydroxy steroid criteria"