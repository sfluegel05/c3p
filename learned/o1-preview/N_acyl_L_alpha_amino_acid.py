"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
"""
Classifies: N-acyl-L-alpha-amino acid
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

    # Parse SMILES string into RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for N-acyl-alpha-amino acid backbone (without chirality constraints)
    pattern = Chem.MolFromSmarts("C(NC(=O)[#6])[#6]C(=O)O")
    if not mol.HasSubstructMatch(pattern):
        return False, "Does not contain N-acyl-alpha-amino acid backbone"

    # Count the number of amide bonds in the molecule
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) > 1:
        return False, f"Contains {len(amide_matches)} amide bonds, may be a peptide"

    return True, "Molecule is an N-acyl-L-alpha-amino acid"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'N-acyl-L-alpha-amino acid',
        'definition': 'Any L-alpha-amino acid carrying an N-acyl substituent.'
    }
}