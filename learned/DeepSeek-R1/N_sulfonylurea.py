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

    # Define urea pattern: two nitrogens connected by a carbonyl
    urea_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[NX3]")
    urea_matches = mol.GetSubstructMatches(urea_pattern)
    if not urea_matches:
        return False, "No urea group detected"

    # Define sulfonyl group attached to nitrogen: N-S(=O)(=O)
    sulfonyl_pattern = Chem.MolFromSmarts("[NX3][SX4](=[OX1])(=[OX1])")
    sulfonyl_matches = mol.GetSubstructMatches(sulfonyl_pattern)
    if not sulfonyl_matches:
        return False, "No sulfonyl group attached to nitrogen"

    # Check if sulfonamide nitrogen is part of any urea group
    for sulf_match in sulfonyl_matches:
        sulf_n_idx = sulf_match[0]  # Nitrogen connected to sulfonyl
        for urea_match in urea_matches:
            # Check if sulfonamide nitrogen is one of the urea nitrogens
            if sulf_n_idx in urea_match[:2]:  # First two atoms in urea pattern are nitrogens
                return True, "Urea group with sulfonyl substitution on nitrogen"

    return False, "Sulfonyl group not attached to urea nitrogen"