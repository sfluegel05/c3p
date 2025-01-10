"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    A 16beta-hydroxy steroid has a hydroxyl group at position 16 of the steroid backbone
    with a beta-configuration.

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
    
    # General pattern for most common steroid backbone
    steroid_pattern = Chem.MolFromSmarts(
        "[#6]1:3(-[C@@H]2CC[C@H]4CC[C@@]3([H])[C@H]4C2)([#6]-[#6]-[#6]-1)-[#6]-[#6]-[#6]")

    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Pattern to match 16beta-hydroxy group
    # Note: Typically, you identify position-specific groups relative to a 'core' identifier on the steroid framework
    beta_hydroxy_16_pattern = Chem.MolFromSmarts(
        "[C@H](O)[CH2]C1CC[C@]2(C)CCC3C[CH2]CC[C@]32C1")

    if not mol.HasSubstructMatch(beta_hydroxy_16_pattern):
        return False, "No 16beta-hydroxy group found"
    
    return True, "Contains a steroid backbone with 16beta-hydroxy group"

# Example usage
print(is_16beta_hydroxy_steroid("C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1C[C@H](O)[C@@H]2O"))