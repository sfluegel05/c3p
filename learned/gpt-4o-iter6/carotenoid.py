"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    A carotenoid is a tetraterpenoid (C40) with characteristic structural features.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flexible carbon count check for carotenoids and derivatives
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (36 <= c_count <= 55):
        return False, f"Carbon count {c_count} not consistent with typical carotenoid structure"

    # Conjugation check: enhanced to detect 3 consecutive C=C bonds allows more flexibility
    conjugated_pattern = Chem.MolFromSmarts("C=C-C=C")
    conjugated_matches = mol.GetSubstructMatches(conjugated_pattern)
    if len(conjugated_matches) < 3:
        return False, f"Conjugation pattern too fragmented; detected {len(conjugated_matches)} sections of extended conjugation"

    # Adaptive detection of cyclization - expanded range to include common five or six-membered rings
    cyclic_pattern_5 = Chem.MolFromSmarts("C1CCCC1")  # Cyclopentane-like
    cyclic_pattern_6 = Chem.MolFromSmarts("C1CCCCC1")  # Cyclohexane-like
    if not (mol.HasSubstructMatch(cyclic_pattern_5) or mol.HasSubstructMatch(cyclic_pattern_6)):
        return False, "No notable cyclization (rings) detected"

    # Check for typical oxygen-containing functional groups (e.g. C-O, C=O)
    any_oxygen_groups = any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())
    # Adaptable to include more diverse carotenoids
    if not any_oxygen_groups:
        return False, "No notable oxygen-based functional groups found"

    return True, "Structure consistent with carotenoid (C40 tetraterpenoid) with extended conjugation and probable modifications"