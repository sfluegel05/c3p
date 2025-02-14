"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    An N-acyl-L-alpha-amino acid is any L-alpha-amino acid carrying an N-acyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for N-acyl-L-alpha-amino acid with L-stereochemistry
    n_acyl_l_alpha_amino_acid = Chem.MolFromSmarts('[C@@H](NC(=O)[#6])[#6]C(=O)O')
    
    # Check for match with L-stereochemistry
    if mol.HasSubstructMatch(n_acyl_l_alpha_amino_acid):
        return True, "Molecule is an N-acyl-L-alpha-amino acid"

    # SMARTS pattern without stereochemistry specification
    n_acyl_alpha_amino_acid = Chem.MolFromSmarts('C(NC(=O)[#6])[#6]C(=O)O')
    
    # Check for match without stereochemistry
    if mol.HasSubstructMatch(n_acyl_alpha_amino_acid):
        return True, "Molecule is an N-acyl-alpha-amino acid (without specified stereochemistry)"

    return False, "Molecule does not contain N-acyl-L-alpha-amino acid"