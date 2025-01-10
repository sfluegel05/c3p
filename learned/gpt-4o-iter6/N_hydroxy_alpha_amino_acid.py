"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
from rdkit import Chem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino-acid based on its SMILES string.
    This requires an amino acid in which at least one hydrogen attached to the amino group 
    is replaced by a hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a N-hydroxy-alpha-amino-acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Alpha amino acid backbone with variance
    amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]-[C;A][C](=O)[O,O-]")  # Flexible primary or secondary nitrogen next to alpha-carbon
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No alpha amino acid backbone found"

    # Search for at least one hydroxy substitution on the nitrogen
    n_hydroxy_variants = "[N;H1]O"  # Single variant
    bidentate_hydroxy = "[NX3](O)O"  # Bidentate variant
    
    n_hydroxy_patterns = [
        Chem.MolFromSmarts(n_hydroxy_variants),
        Chem.MolFromSmarts(bidentate_hydroxy)
    ]

    if not any(mol.HasSubstructMatch(pattern) for pattern in n_hydroxy_patterns):
        return False, "No N-hydroxy modification found"
    
    return True, "Contains amino acid backbone with N-hydroxy modification"

__metadata__ = {  
    'chemical_class': {   
        'name': 'N-hydroxy-alpha-amino-acid',
        'definition': 'Any amino acid in which at least one hydrogen attached to the amino group is replaced by a hydroxy group.',
    },
    'message': None,
    'attempt': 2,
    'success': False,
}