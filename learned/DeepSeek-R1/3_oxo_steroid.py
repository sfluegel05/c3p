"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: 3-oxo steroid (CHEBI:37781)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid has a steroid nucleus with an oxo group at position 3.

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

    # Define steroid nucleus pattern (four fused rings: three 6-membered, one 5-membered)
    steroid_smarts = Chem.MolFromSmarts("[*]1@[*]@[*]2@[*]@[*]3@[*]@[*]4@[*]@[*]5@[*]@[*]1@[*]@[*]@[*]@[*]2@[*]@[*]@[*]3@[*]@[*]@[*]4@5")
    if not mol.HasSubstructMatch(steroid_smarts):
        return False, "No steroid nucleus detected"

    # Check for oxo group at position 3 (specific position in the A-ring)
    # Position 3 is typically part of the first 6-membered ring (A-ring)
    oxo_pattern = Chem.MolFromSmarts("[C]=[O]")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if not oxo_matches:
        return False, "No oxo group found"

    # Verify at least one oxo group is at position 3
    # This requires knowing the steroid atom numbering, which is complex
    # Alternative approach: Check if any carbonyl is in the A-ring vicinity
    a_ring_smarts = Chem.MolFromSmarts("[*]1@[*]@[*]@[*]@[*]@[*]1")  # A-ring (6-membered)
    a_ring_matches = mol.GetSubstructMatches(a_ring_smarts)
    for oxo_atom in oxo_matches:
        for a_ring in a_ring_matches:
            if oxo_atom[0] in a_ring:
                # Assuming position 3 is the third atom in the A-ring SMARTS match
                if a_ring.index(oxo_atom[0]) == 2:
                    return True, "3-oxo group found on steroid nucleus"

    return False, "Oxo group not found at position 3"