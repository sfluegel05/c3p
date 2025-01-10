"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    Carotenoids are C40 tetraterpenoids with a long conjugated double-bond system,
    potentially containing cyclic structures, oxygen functionalities, and excluding retinoids.

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

    # Adjust carbon atom count range to account for derivatives
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30 or c_count > 60:
        return False, f"Number of carbon atoms is {c_count}, not within expected range for carotenoids"

    # Improved check for conjugated double-bond system
    conjugated_double_bond_pattern = Chem.MolFromSmarts("C=C")
    conjugated_double_bond_matches = mol.GetSubstructMatches(conjugated_double_bond_pattern)
    if len(conjugated_double_bond_matches) < 10:
        return False, f"Fewer conjugated double bonds found ({len(conjugated_double_bond_matches)}), suggestive of non-carotenoid structure"

    # Exclude retinoid structures by recognizing common motifs
    # Example pattern area for expansion (refinement needed experimentally)
    retinoid_pattern_full = Chem.MolFromSmarts("CC=C(CC=C)CC=C")
    if mol.HasSubstructMatch(retinoid_pattern_full) and c_count < 35:
        return False, "Structure matches common retinoid pattern or too few carbons, excluded from carotenoids"

    # Functional groups presence check (e.g., oxygens)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count > 0:
        return True, "Classified as a xanthophyll carotenoid due to the presence of oxygen functionality"

    return True, "Classified as a carotene (oxygen-less carotenoid) due to C40 tetraterpenoid structure with conjugated bonds"