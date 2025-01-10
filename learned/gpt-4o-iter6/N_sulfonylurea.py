"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
from rdkit import Chem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    An N-sulfonylurea contains a urea group where one nitrogen is substituted 
    with a sulfonyl group (-S(=O)(=O)-).

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
    # N-sulfonyl: N-S(=O)(=O)
    nsulfonyl_pattern = Chem.MolFromSmarts("N[SX4](=O)(=O)")

    # Find matches of the urea pattern
    urea_matches = mol.GetSubstructMatches(urea_pattern)
    if not urea_matches:
        return False, "No urea group found"
    
    # Find matches of the sulfonyl pattern
    nsulfonyl_matches = mol.GetSubstructMatches(nsulfonyl_pattern)
    if not nsulfonyl_matches:
        return False, "No N-sulfonyl group found"

    # We need to verify that the sulfonyl is attached to one of the urea nitrogen atoms
    for urea_match in urea_matches:
        # The nitrogen atoms are at indices 0 and 2 of urea_match
        urea_nitrogens = [urea_match[0], urea_match[2]]
        for n_idx in urea_nitrogens:
            nitrogen = mol.GetAtomWithIdx(n_idx)
            for neighbor in nitrogen.GetNeighbors():
                if neighbor.GetIdx() in [match[1] for match in nsulfonyl_matches]:
                    return True, "Contains N-sulfonylurea moiety"
    
    return False, "No N-sulfonyl substitution on urea nitrogen"