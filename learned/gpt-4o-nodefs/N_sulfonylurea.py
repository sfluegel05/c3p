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
    
    # Improved pattern to match sulfonyl group attached to nitrogen
    sulfonyl_nitrogen_pattern = Chem.MolFromSmarts("S(=O)(=O)N")
    if not mol.HasSubstructMatch(sulfonyl_nitrogen_pattern):
        return False, "No N-sulfonyl group found"
    
    # Improved pattern to match broader urea moiety configurations
    urea_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[NX3]")
    if not mol.HasSubstructMatch(urea_pattern):
        return False, "No urea moiety found"

    return True, "Contains both N-sulfonyl group and urea moiety"