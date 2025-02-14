"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1)-steroid based on its SMILES string.
    A 3-oxo-Delta(1)-steroid is a steroid with a ketone at position 3 and a double bond
    between positions 1 and 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(1)-steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a single SMARTS pattern for the steroid core, 3-oxo group and delta(1) double bond
    # This ensures all conditions are met in the same ring system.
    # The numbering refers to specific positions in the steroid core
    combined_pattern = Chem.MolFromSmarts("[C]1=[C][C]2[C]3[C]([C](=[O])1)([C][C]4[C]3[C]([C]2)[C]4)")
    
    if not mol.HasSubstructMatch(combined_pattern):
        return False, "Molecule does not match the combined 3-oxo-Delta(1)-steroid pattern."

    return True, "Molecule is a 3-oxo-Delta(1)-steroid."