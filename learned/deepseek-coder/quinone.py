"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    A quinone is a compound with a fully conjugated cyclic dione structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a cyclic structure
    if not mol.GetRingInfo().NumRings():
        return False, "No cyclic structure found"

    # Look for two carbonyl groups (C=O) in the molecule
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if len(carbonyl_matches) < 2:
        return False, f"Found {len(carbonyl_matches)} carbonyl groups, need at least 2"

    # Check if the carbonyl groups are in the same ring system
    # Create a pattern that looks for two carbonyls in any ring system
    ring_system_pattern = Chem.MolFromSmarts("[CX3]=[OX1]~*~*~[CX3]=[OX1]")
    if not mol.HasSubstructMatch(ring_system_pattern):
        return False, "No two carbonyl groups in the same ring system"

    # Check that the carbonyl groups are in a conjugated system
    # Use a more flexible pattern to detect conjugation
    conjugated_pattern = Chem.MolFromSmarts("[CX3]=[OX1]~[CX3]=,*~[CX3]=[OX1]")
    if not mol.HasSubstructMatch(conjugated_pattern):
        return False, "No conjugated system between carbonyl groups"

    return True, "Contains a fully conjugated cyclic dione structure"