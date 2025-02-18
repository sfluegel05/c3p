"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
"""
Classifies: N-sulfonylurea
"""
from rdkit import Chem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    An N-sulfonylurea is a urea where one nitrogen is substituted with a sulfonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-sulfonylurea, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Combined pattern for N-sulfonylurea core: urea with sulfonyl directly attached to one nitrogen
    # Matches either N-SO2-C(=O)-N or N-C(=O)-N-SO2 arrangements
    pattern = "C(=O)(N[SD4](=O)(=O))N"
    core_pattern = Chem.MolFromSmarts(pattern)
    if core_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check for presence of the core structure
    has_core = mol.HasSubstructMatch(core_pattern)
    if has_core:
        return True, "Contains urea group with sulfonyl substitution on nitrogen"

    # Alternative check using expanded SMARTS to capture all possible arrangements
    # This covers both possible positions for the sulfonyl group relative to urea
    alt_pattern = Chem.MolFromSmarts("[NX3]C(=O)[NX3][SD4](=O)(=O)")
    if alt_pattern and mol.HasSubstructMatch(alt_pattern):
        return True, "Contains urea group with sulfonyl substitution on nitrogen"

    # Additional check for sulfonyl group attached to the other urea nitrogen
    alt_pattern2 = Chem.MolFromSmarts("[NX3][SD4](=O)(=O)C(=O)[NX3]")
    if alt_pattern2 and mol.HasSubstructMatch(alt_pattern2):
        return True, "Contains urea group with sulfonyl substitution on nitrogen"

    # Verify presence of urea group (as fallback)
    urea_pattern = Chem.MolFromSmarts("[NX3]C(=O)[NX3]")
    if not mol.HasSubstructMatch(urea_pattern):
        return False, "No urea group detected"

    # Verify sulfonamide group connected to urea nitrogen
    sulfonyl_n_pattern = Chem.MolFromSmarts("[NX3][SD4](=O)(=O)")
    sulfonyl_matches = mol.GetSubstructMatches(sulfonyl_n_pattern)
    urea_matches = mol.GetSubstructMatches(urea_pattern)

    # Check if any sulfonamide nitrogen is part of the urea group
    urea_n_indices = set()
    for match in urea_matches:
        urea_n_indices.update([match[0], match[2]])  # Positions of N atoms in urea

    for sulf_match in sulfonyl_matches:
        sulf_n_idx = sulf_match[0]
        if sulf_n_idx in urea_n_indices:
            return True, "Urea group with sulfonyl substitution on nitrogen"

    return False, "Sulfonyl group not attached to urea nitrogen"