"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
from rdkit import Chem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    A 17beta-hydroxy steroid must have a hydroxy group at position 17 in beta-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flexible steroid backbone pattern
    steroid_patterns = [
        Chem.MolFromSmarts('[#6]1[#6][#6][#6]2[#6][#6][#6]3[#6][#6][#6][#6]4[#6][#6][#6][#6]4[#6]3[#6]2[#6]1'),  # Simplified steroid backbone
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in steroid_patterns):
        return False, "No steroid backbone found"

    # Pattern to identify 17beta-hydroxy group
    hydroxy_17beta_pattern = Chem.MolFromSmarts('[C@@H]1([C@@](O)(C2)C3CC[C@@H]2C4CC[C@@H](C4)C3)[#6]C([#6])[#6]1')
    if not mol.HasSubstructMatch(hydroxy_17beta_pattern):
        return False, "No 17beta-hydroxy group in beta-configuration found"
    
    return True, "Contains 17beta-hydroxy group with steroid backbone configuration"