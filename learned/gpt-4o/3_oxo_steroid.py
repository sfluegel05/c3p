"""
Classifies: CHEBI:47788 3-oxo steroid
"""
from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is characterized by an oxo (ketone) group at the third position on the steroid backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general SMARTS pattern for steroid backbone (rings: B-C-D)
    steroid_backbone_pattern = Chem.MolFromSmarts("[#6]1[#6][#6]2[#6][#6]3[#6][#6][#6]4[#6][#6][#6]1[#6]2[#6][#6]3")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No steroid backbone found"

    # Identify the correct 3-position for the oxo group in the steroid structure
    # Use a pattern that identifies a C=O group at possible positions
    oxo_group_pattern = Chem.MolFromSmarts("[#6;r3](=O)")
    if not any(mol.GetSubstructMatches(oxo_group_pattern)):
        return False, "No 3-oxo group found"

    return True, "Molecule is a 3-oxo steroid"