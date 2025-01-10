"""
Classifies: CHEBI:27300 vitamin D
"""
"""
Classifies: vitamin D
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    Vitamin D compounds are secosteroids with a broken B-ring, resulting in three rings.
    They have a conjugated triene system and contain hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a vitamin D compound, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings != 3:
        return False, f"Molecule has {num_rings} rings, expected 3 for a secosteroid"

    # Check for conjugated triene system (C=C-C=C-C=C)
    triene_pattern = Chem.MolFromSmarts('C=C-C=C-C=C')
    if not mol.HasSubstructMatch(triene_pattern):
        return False, "No conjugated triene system found"

    # Check for at least one hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    hydroxyls = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyls) == 0:
        return False, "No hydroxyl groups found"

    # Optionally, check for secosteroid core using a subset of the steroid skeleton
    # This pattern represents a fused ring system with one ring opened
    secosteroid_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CC=C2')
    if not mol.HasSubstructMatch(secosteroid_pattern):
        return False, "Secosteroid core structure not found"

    # Additional checks can be added as needed to improve specificity

    return True, "Molecule matches vitamin D structural features"