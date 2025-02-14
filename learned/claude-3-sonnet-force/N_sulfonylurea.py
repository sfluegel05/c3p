"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
"""
Classifies: CHEBI:50123 N-sulfonylurea
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    An N-sulfonylurea is a urea in which one of the hydrogens attached to a nitrogen
    of the urea group is replaced by a sulfonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-sulfonylurea, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for urea group (-N(C=O)N-)
    urea_pattern = Chem.MolFromSmarts("[NH1]C(=O)N")
    urea_matches = mol.GetSubstructMatches(urea_pattern)
    if not urea_matches:
        return False, "No urea group found"

    # Look for sulfonyl group (-SO2-) attached to a nitrogen of the urea
    sulfonylurea_pattern = Chem.MolFromSmarts("N(S(=O)(=O))C(=O)N")
    sulfonylurea_matches = mol.GetSubstructMatches(sulfonylurea_pattern)
    if not sulfonylurea_matches:
        return False, "No N-sulfonylurea group found"

    # Check if the sulfonylurea group is part of the urea group
    for urea_match, sulfonylurea_match in zip(urea_matches, sulfonylurea_matches):
        if set(urea_match) & set(sulfonylurea_match):
            return True, "Contains an N-sulfonylurea group (urea with a sulfonyl group attached to a nitrogen)"

    return False, "N-sulfonylurea group not part of urea group"