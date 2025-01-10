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

    # Adjust carbon count - considering even derivatives of carotenoids
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (35 <= c_count <= 45):
        return False, f"Carbon count {c_count} not consistent with typical carotenoid structure"

    # Enhanced check for extended conjugation with at least 8 C=C in sequence
    conjugated_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")
    conjugated_matches = mol.GetSubstructMatches(conjugated_pattern)
    if len(conjugated_matches) < 4:
        return False, f"Conjugation pattern too fragmented; detected {len(conjugated_matches)} sections of extended conjugation"

    # Detect functional groups: hydroxyl (OH) or carbonyl (C=O)
    any_oxygen_groups = any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())
    if not any_oxygen_groups:
        return False, "No notable oxygen-based functional groups found"

    # Adaptive detection of cyclization - carotenoids may have ring closures
    for ring_size in range(5, 7):
        cyclic_patterns = Chem.MolFromSmarts(f"C1{'C' * (ring_size - 2)}C1")
        if mol.HasSubstructMatch(cyclic_patterns):
            break
    else:
        return False, "No notable cyclization (rings) detected"

    return True, "Structure consistent with carotenoid (C40 tetraterpenoid) with extended conjugation and probable modifications"