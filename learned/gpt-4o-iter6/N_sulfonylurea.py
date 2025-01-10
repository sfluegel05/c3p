"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
from rdkit import Chem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    A true N-sulfonylurea contains a urea group with one nitrogen atom replaced
    by a sulfonyl group.

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
    # N-sulfonyl pattern (NSO2): [N;X3][SX4](=O)(=O) (only attached nitrogen)
    nsulfonyl_pattern = Chem.MolFromSmarts("[NX3;!$(NC=O)][SX4](=O)(=O)")

    # Find matches of the urea pattern
    urea_matches = mol.GetSubstructMatches(urea_pattern)
    if not urea_matches:
        return False, "No urea group found"

    # Find all matches of the correct N-sulfonyl pattern
    nsulfonyl_matches = mol.GetSubstructMatches(nsulfonyl_pattern)
    nsulfonyl_n_index = [match[0] for match in nsulfonyl_matches]

    # Check if any urea nitrogen is part of the correct N-sulfonyl group
    for urea_match in urea_matches:
        # Nitrogen atoms are at indices 0 and 2 of urea_match
        if urea_match[0] in nsulfonyl_n_index or urea_match[2] in nsulfonyl_n_index:
            return True, "Contains N-sulfonylurea moiety"
    
    return False, "No N-sulfonyl substitution on urea nitrogen"