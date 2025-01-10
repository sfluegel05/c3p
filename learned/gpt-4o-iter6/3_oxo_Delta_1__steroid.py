"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_1_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid based on its SMILES string.
    A 3-oxo-Delta(1) steroid is characterized by a steroid structure with
    a double bond between positions 1 and 2 and a ketone group at the 3-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(1) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for a ketone group at the 3-position
    ketone_pattern = Chem.MolFromSmarts("[C;R1](=O)[C;R1]")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group at the 3-position"

    # SMARTS pattern for a double bond between positions 1 and 2
    delta_1_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(delta_1_bond_pattern):
        return False, "No double bond between positions 1 and 2"
    
    # SMARTS pattern for a general steroid structure (ABCDE ring system)
    # This is more flexible and may need further adjustment based on stereo descriptors
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CC2CCC3C4C(CCC4)CCC3C2C1")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No valid steroid backbone found"
    
    return True, "Matches the 3-oxo-Delta(1) steroid pattern"