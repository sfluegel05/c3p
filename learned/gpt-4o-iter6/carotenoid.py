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

    # Carotenoids are derived from a C40 framework, flexible around 40 carbons
    # Tetraterpenoid typically implies ~40 carbons or conmplex derivative structures
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (33 <= c_count <= 48):  # Flexible range based on derivations
        return False, f"Carbon count {c_count} not typical for carotenoid derivatives"

    # Extended conjugated system: at least several consecutive conjugated C=C bonds
    conjugated_pattern = Chem.MolFromSmarts("C=C(-C=C)+")
    conjugated_matches = mol.GetSubstructMatches(conjugated_pattern)
    if len(conjugated_matches) < 5:
        return False, f"Extended conjugation pattern lacking; detected {len(conjugated_matches)} sections"

    # Looser cyclization check for more flexibility - explore multiple ring types or lack thereof
    cyclization_flexible_patterns = [
        Chem.MolFromSmarts("C1CCCC1"),  # Cyclopentane-like
        Chem.MolFromSmarts("C1CCCCC1"),  # Cyclohexane-like
        Chem.MolFromSmarts("C1C=CCC1"),  # Cyclopentene-like
        Chem.MolFromSmarts("C1C=CCCC1")  # Cyclohexene-like
    ]
    # Evaluate if any acceptable cyclization is present or can be derived
    has_ring = any(mol.HasSubstructMatch(pattern) for pattern in cyclization_flexible_patterns)
    if not has_ring and c_count == 40:  # Allow more flexibility if precisely a C40 backbone
        return True, "Structure consistent with carotenoid (C40 tetraterpenoid) despite limited cyclization"

    # Allow presence of oxygen-based functional groups (reflects class diversity by including xanthophylls)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count > 0:
        return True, "Presence of oxygen groups indicates possible carotenoid derivative (xanthophyll or similar)"

    # Final fallback return if comprehensive traits match
    return True, "Structure generally consistent with carotenoid framework after pattern analysis"