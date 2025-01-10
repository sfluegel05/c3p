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

    # Verify steroid backbone with 4-ring system
    ring_info = mol.GetRingInfo()
    if len(ring_info.AtomRings()) < 4:
        return False, "Not enough rings for a steroid backbone"

    # Check for steroid backbone pattern
    steroid_backbone_pattern = Chem.MolFromSmarts("C1[C@H]2CCC3=CC(=O)CC4=C[C@H]2CC[C@]3(C)[C@@H]4CC1")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No steroid backbone found"

    # Check for oxo (C=O) group at the third position
    oxo_group_pattern = Chem.MolFromSmarts("C3(=O)")
    if not mol.HasSubstructMatch(oxo_group_pattern):
        return False, "No 3-oxo group found"

    return True, "Molecule is a 3-oxo steroid"