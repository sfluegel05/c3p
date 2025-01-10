"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    This class is defined by an L-alpha-amino acid backbone with an N-acyl substituent.

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

    # Look for L-alpha-amino acid backbone pattern: chiral center with amino group and carboxylic acid
    amino_acid_pattern = Chem.MolFromSmarts("[C@H](N)C(=O)O")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No L-alpha-amino acid backbone found"

    # Look for N-acyl group: N attached to C=O
    n_acyl_pattern = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(n_acyl_pattern):
        return False, "No N-acyl group found"

    return True, "Contains L-alpha-amino acid backbone with N-acyl substituent"