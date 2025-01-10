"""
Classifies: CHEBI:76983 N-sulfonylurea
"""
from rdkit import Chem

def is_N_sulfonylurea(smiles: str):
    """
    Determines if a molecule is an N-sulfonylurea based on its SMILES string.
    An N-sulfonylurea contains a sulfonyl group (-SO2-) attached to a nitrogen atom,
    and a urea moiety (-NH-C(=O)-NH-) in its structure, specifically structured.
    
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
    
    # Pattern for N-sulfonylurea: Ensure N-sulfonyl is attached directly to urea nitrogen
    # Contains the pattern: [NX3]S(=O)(=O)N[C](=O)N or an adjacent variant
    # The urea nitrogen is part of the sulfonyl group link
    n_sulfonylurea_pattern = Chem.MolFromSmarts("N[S](=O)(=O)NC(=O)N")
    
    if mol.HasSubstructMatch(n_sulfonylurea_pattern):
        return True, "Contains N-sulfonyl group directly bonded to urea moiety"
    
    return False, "Missing specific N-sulfonylurea linkage in the structure"