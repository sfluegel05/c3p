"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: N-hydroxy-alpha-amino-acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino acid based on its SMILES string.
    These compounds have a carboxylic acid group and at least one hydrogen attached 
    to the amino group is replaced by a hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-hydroxy-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Define patterns for N-hydroxy groups
    n_hydroxy_patterns = [
        # Simple N-hydroxy
        "[C][N]O",
        # N,N-dihydroxy
        "[C][N](O)O",
        # Hydroxyimino groups (various configurations)
        "[C]N=NO",
        "[C]NC(=NO)",
        "[C]NC(=N)NO",
        # Cover both E and Z isomers
        "[C]N/C(=N/O)",
        "[C]N\C(=N\O)",
        "[C]N/C(=N\O)",
        "[C]N\C(=N/O)",
        # Additional patterns for completeness
        "[C]N([H])O",
        "[C]N(O)[H]"
    ]
    
    # Check for alpha-amino acid backbone with N-hydroxy group
    found_n_hydroxy = False
    matched_pattern = None
    
    for pattern in n_hydroxy_patterns:
        substructure = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(substructure):
            found_n_hydroxy = True
            matched_pattern = pattern
            break
            
    if not found_n_hydroxy:
        return False, "No N-hydroxy group found"

    # Check for alpha-amino acid backbone
    # More general pattern that captures various configurations
    alpha_patterns = [
        # Basic alpha-amino acid pattern
        "[C]([C](=O)O)([#1,*])[N,n]",
        # Pattern for cyclic amino acids
        "C1([C](=O)O)([#1,*])N([OH1,OH0])CCC1",
        # Pattern for substituted alpha carbons
        "[C]([C](=O)O)([#1,*])[N]"
    ]
    
    found_alpha = False
    for pattern in alpha_patterns:
        substructure = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(substructure):
            found_alpha = True
            break
    
    if not found_alpha:
        return False, "Not an alpha-amino acid structure"

    # Additional validation for correct connectivity
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    n_hydroxy_matches = mol.GetSubstructMatches(Chem.MolFromSmarts(matched_pattern))
    
    if not carboxyl_matches or not n_hydroxy_matches:
        return False, "Required groups not properly connected"

    # Success cases with specific classifications
    if "[N](O)O" in matched_pattern:
        return True, "Contains N,N-dihydroxy group on alpha-amino acid"
    elif "N=NO" in matched_pattern or "C(=NO)" in matched_pattern:
        return True, "Contains hydroxyimino group on alpha-amino acid"
    else:
        return True, "Contains N-hydroxy group on alpha-amino acid"