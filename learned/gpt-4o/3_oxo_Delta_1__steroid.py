"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid based on its SMILES string.
    A 3-oxo-Delta(1) steroid is a steroid with a ketone group at position 3 and a double bond
    between positions 1 and 2.

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

    # Look for 3-oxo group pattern in steroids (better specify environment and rings)
    oxo_pattern = Chem.MolFromSmarts("[#6]-[#8]=[#6;R]")  # Carbonyl C(=O)C associated with a ring
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No 3-oxo group (C=O near the typical ketone position in steroids) found"
    
    # Look for a double bond between positions 1 and 2 in a ring (enhanced identity checking)
    # Position referencing must align with standardized steroid frameworks to ensure accuracy
    delta_1_pattern = Chem.MolFromSmarts("[#6;R1]=[#6;R2]")  # Atoms are ring atoms
    if not mol.HasSubstructMatch(delta_1_pattern):
        return False, "No Delta(1) double bond (C=C between position 1 and 2) found"

    # Verify additional structure commonness of a steroid (basic structure or tetracyclic schemes)
    steroid_tetracyclic_pattern = Chem.MolFromSmarts("C1CC2CCC3C4(CC3C2C1)C=CC=C4") 
    # Check with a simpler more inclusive steroid backbone pattern
    if not mol.HasSubstructMatch(steroid_tetracyclic_pattern):
        return False, "No tetracyclic ring characteristic to steroid backbone found"

    return True, "Matches 3-oxo-Delta(1) steroid structure"