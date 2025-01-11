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

    # Define a comprehensive SMARTS pattern for N-sulfonylurea
    nsulfonylurea_pattern = Chem.MolFromSmarts("NC(=O)N[SX4](=O)(=O)")

    # Find matches of the N-sulfonylurea pattern
    nsulfonylurea_matches = mol.GetSubstructMatches(nsulfonylurea_pattern)
    if nsulfonylurea_matches:
        return True, "Contains N-sulfonylurea moiety"

    return False, "No N-sulfonyl substitution on urea nitrogen"