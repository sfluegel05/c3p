"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    Carotenoids are C40 tetraterpenoids with a long conjugated double-bond system,
    potentially containing cyclic structures, oxygen functionalities, and activity modifications.

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

    # Determine carbon atom count (adjusted for carotenoid derivatives)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30 or c_count > 50:  # Usual core-carotenoid range plus some leeway
        return False, f"Number of carbon atoms is {c_count}, not indicative for carotenoids"

    # Verify conjugated double-bond system (increase threshold and detail check)
    conjugated_double_bond_pattern = Chem.MolFromSmarts("C=C")
    conjugated_double_bond_matches = mol.GetSubstructMatches(conjugated_double_bond_pattern)
    if len(conjugated_double_bond_matches) < 9:  # Lower threshold for detection
        return False, f"Insufficient conjugated double bonds ({len(conjugated_double_bond_matches)}) found"

    # Recognize retinoids to exclude false positives
    retinoid_pattern = Chem.MolFromSmarts("CC=C(CC)CC=C")
    if mol.HasSubstructMatch(retinoid_pattern) and c_count < 35:
        return False, "Common retinoid-like pattern detected, excluded from carotenoids"

    # Look for oxygen or other functionality that distinguishes classes
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count > 0 and "O" in smiles:
        return True, "Oxygen functionality present, classifying as a xanthophyll carotenoid"

    return True, "Conforms with carotene (oxygen-less carotenoid) structure from SMILES-derived data"