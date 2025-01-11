"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
from rdkit import Chem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    An N-sulfonylurea contains a sulfonyl group (-SO2-) attached to a nitrogen atom,
    and a urea moiety (-NH-C(=O)-NH-) in its structure.

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

    # Look for sulfonyl group (S(=O)(=O)) attached to a nitrogen
    sulfonyl_nitrogen_pattern = Chem.MolFromSmarts("S(=O)(=O)[N]")
    if not mol.HasSubstructMatch(sulfonyl_nitrogen_pattern):
        return False, "No N-sulfonyl group found"
    
    # Look for urea moiety (-NH-C(=O)-NH-)
    urea_pattern = Chem.MolFromSmarts("N-C(=O)-N")
    if not mol.HasSubstructMatch(urea_pattern):
        return False, "No urea moiety found"

    return True, "Contains both N-sulfonyl group and urea moiety"