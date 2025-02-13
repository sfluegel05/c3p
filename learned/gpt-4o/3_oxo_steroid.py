"""
Classifies: CHEBI:47788 3-oxo steroid
"""
from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is characterized by a steroidal framework with an oxo group (C=O) at the 3rd position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid backbone pattern
    # (4 fused rings with allowance for common variations)
    steroid_pattern = Chem.MolFromSmarts('[#6]1([#6])[#6][#6]2[#6]3[#6][#8]C(=O)[#6]=[C]4[#6]3CC[C@]12C[C@@H](O)[C@]4(C)C')
    
    # Check for steroid backbone
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Define the oxo group pattern at position 3
    # A more refined, specific pattern is needed here
    # Identify the oxo group in context of known steroid structure
    oxo_pattern = Chem.MolFromSmarts('[#8]=[C]3CC[C@](C)(O)[#6]C(=O)[#6]4[C@@H]3C[C@]12CC')
    
    # Check if this specific pattern for a 3-oxo structure matches
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No oxo group (C=O) at the 3rd position found"

    return True, "Contains a steroid backbone with a 3-oxo group at the correct position."