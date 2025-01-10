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

    # Check for L-alpha-amino acid backbone: Include central carbon with attached NH2 and C(=O)OH
    l_alpha_amino_acid_pattern = Chem.MolFromSmarts("[C@@H](N)(*)C(=O)O")
    if not mol.HasSubstructMatch(l_alpha_amino_acid_pattern):
        return False, "No L-alpha-amino acid backbone found with correct stereochemistry"

    # Check for N-acyl group: Use flexible acyl pattern wide-encompassing various R-groups
    n_acyl_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
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