"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
from rdkit import Chem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is a N-hydroxy-alpha-amino-acid based on its SMILES string.
    This requires an amino acid in which at least one hydrogen attached to the amino group 
    is replaced by a hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a N-hydroxy-alpha-amino-acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for the alpha amino acid backbone: [C](C(=O)O)[N]
    amino_acid_pattern = Chem.MolFromSmarts("C[C@H](N)C(=O)O")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No alpha amino acid backbone found"
    
    # Look for N-hydroxy modification: [N](-[O])(|[O])
    n_hydroxy_pattern = Chem.MolFromSmarts("[NX3](-[OX1])[OX1]")
    n_hydroxy_matches = mol.GetSubstructMatches(n_hydroxy_pattern)
    if len(n_hydroxy_matches) == 0:
        return False, "No N-hydroxy modification found"

    return True, "Contains amino acid backbone with N-hydroxy modification"

__metadata__ = {  
    'chemical_class': {   
        'name': 'N-hydroxy-alpha-amino-acid',
        'definition': 'Any amino acid in which at least one hydrogen attached to the amino group is replaced by a hydroxy group.',
    },
    'message': None,
    'attempt': 0,
    'success': True,
}