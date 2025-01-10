"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    An N-acyl-L-alpha-amino acid is defined as an L-alpha-amino acid carrying an N-acyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for L-alpha-amino acid backbone
    # Pattern captures a C-alpha center bound to an amine and a carboxyl group, with a specific stereochemistry
    l_alpha_amino_acid_pattern = Chem.MolFromSmarts("[C@H](N)(*)C(=O)O")
    if not mol.HasSubstructMatch(l_alpha_amino_acid_pattern):
        return False, "No L-alpha-amino acid backbone with correct stereochemistry found"

    # Check for N-acyl group
    # Flexible acyl pattern with nitrogen carrying the substituent
    n_acyl_pattern = Chem.MolFromSmarts("N[C](=O)C")
    if not mol.HasSubstructMatch(n_acyl_pattern):
        return False, "No N-acyl group found attached to nitrogen"

    # If all patterns matched
    return True, "Contains a valid L-alpha-amino acid backbone with an N-acyl substituent"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:17855',
        'name': 'N-acyl-L-alpha-amino acid',
        'definition': 'Any L-alpha-amino acid carrying an N-acyl substituent.'
    }
}