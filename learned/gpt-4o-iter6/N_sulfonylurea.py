"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
from rdkit import Chem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    A true N-sulfonylurea contains a urea moiety with one nitrogen atom
    substituted by a sulfonyl group (-S(=O)(=O)-).

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

    # Define SMARTS patterns
    # Urea group: N-C(=O)-N
    urea_pattern = Chem.MolFromSmarts("NC(=O)N")
    # More general pattern for N-sulfonyl group: N-S(=O)(=O)
    nsulfonyl_pattern = Chem.MolFromSmarts("N[SX4](=O)(=O)")

    # Find matches of the urea pattern
    urea_matches = mol.GetSubstructMatches(urea_pattern)
    if not urea_matches:
        return False, "No urea group found"

    # Find all matches of the N-sulfonyl pattern
    nsulfonyl_matches = mol.GetSubstructMatches(nsulfonyl_pattern)
    if not nsulfonyl_matches:
        return False, "No N-sulfonyl group found"

    # Map indexes for N-sulfonyl nitrogen
    nsulfonyl_n_index = [match[0] for match in nsulfonyl_matches]

    # Check if any urea nitrogen is part of an N-sulfonyl group
    for urea_match in urea_matches:
        urea_n_atoms = [urea_match[0], urea_match[2]]  # Nitrogen atoms in urea
        if any(n in nsulfonyl_n_index for n in urea_n_atoms):
            return True, "Contains N-sulfonylurea moiety"

    return False, "No N-sulfonyl substitution on urea nitrogen"