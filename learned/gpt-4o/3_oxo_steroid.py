"""
Classifies: CHEBI:47788 3-oxo steroid
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

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

    # Define SMARTS patterns
    steroid_backbone_pattern = Chem.MolFromSmarts("[#6]1[#6][#6]2[#6]([#6]1)[#6][#6]3[#6]([#6]2)[#6][#6]([#6](=O)[#6]3)")
    
    # Verify steroid backbone with at least 3 rings
    ring_info = mol.GetRingInfo()
    if not ring_info.IsInitialized() or len(ring_info.AtomRings()) < 3:
        return False, "Not enough rings for a steroid backbone"

    # Check for matching steroid backbone
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No steroid backbone with 3-oxo group at position 3 found"

    # Confirm 3-oxo group specifically
    # Search for C=O at position 3 of a steroid backbone
    oxo_group_pattern = Chem.MolFromSmarts("C=O")
    if not any(mol.HasSubstructMatch(oxo_group_pattern) for _ in ring_info.AtomRings()):
        return False, "No 3-oxo group at appropriate position"

    return True, "Molecule is a 3-oxo steroid"