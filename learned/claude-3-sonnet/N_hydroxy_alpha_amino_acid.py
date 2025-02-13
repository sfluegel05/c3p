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
    These compounds have a carboxylic acid group and either an N-hydroxy or N,N-dihydroxy group
    attached to an alpha carbon.

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
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Match patterns for N-hydroxy and N,N-dihydroxy groups attached to carbon
    n_hydroxy_pattern = Chem.MolFromSmarts("[CX4][NX3H0]([OX2H1])[H]")  # -C-N(O)H
    n_dihydroxy_pattern = Chem.MolFromSmarts("[CX4][NX3H0]([OX2H1])[OX2H1]")  # -C-N(O)O
    
    has_n_hydroxy = mol.HasSubstructMatch(n_hydroxy_pattern)
    has_n_dihydroxy = mol.HasSubstructMatch(n_dihydroxy_pattern)
    
    if not (has_n_hydroxy or has_n_dihydroxy):
        return False, "No N-hydroxy or N,N-dihydroxy group found"
    
    # Verify the N-hydroxy/dihydroxy group is attached to the alpha carbon
    # First find all carboxyl carbons
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    # For each carboxyl carbon, check if it's connected to a carbon with N-hydroxy group
    valid_alpha = False
    for match in carboxyl_matches:
        carboxyl_carbon = mol.GetAtomWithIdx(match[0])
        for neighbor in carboxyl_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() != 6:  # not carbon
                continue
            # Check if this carbon has N-hydroxy or N,N-dihydroxy group
            if mol.HasSubstructMatch(n_hydroxy_pattern, {neighbor.GetIdx()}):
                valid_alpha = True
                break
            if mol.HasSubstructMatch(n_dihydroxy_pattern, {neighbor.GetIdx()}):
                valid_alpha = True
                break
    
    if not valid_alpha:
        return False, "N-hydroxy/dihydroxy group not attached to alpha carbon"

    hydroxy_type = "N,N-dihydroxy" if has_n_dihydroxy else "N-hydroxy"
    return True, f"Contains {hydroxy_type} group attached to alpha carbon of amino acid"