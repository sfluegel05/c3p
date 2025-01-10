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

    # Define a SMARTS pattern for an alpha amino acid backbone
    # This includes checking for C alpha position with a carboxylic acid and potential substitutions at the nitrogen
    amino_acid_pattern = Chem.MolFromSmarts("N[C@@H](C)C(=O)O")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No alpha amino acid backbone found"

    # Define SMARTS for N-hydroxy pattern, ensuring -N(OH) linkage
    n_hydroxy_pattern = Chem.MolFromSmarts("[N;!R](O)")
    if not mol.HasSubstructMatch(n_hydroxy_pattern):
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