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

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for N-acyl-L-alpha-amino acid
    # This pattern looks for:
    # - An alpha carbon with L-configuration ([C@@H])
    # - Attached to a carboxylic acid group (C(=O)O)
    # - Attached to an N-acylated nitrogen (N-C(=O)-C)
    pattern = Chem.MolFromSmarts('[C@@H](NC(=O)C)[CH2,C](C(=O)O)')

    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Molecule is an N-acyl-L-alpha-amino acid"

    # If not found, check without stereochemistry
    pattern_no_stereo = Chem.MolFromSmarts('C(NC(=O)C)[CH2,C](C(=O)O)')

    if mol.HasSubstructMatch(pattern_no_stereo):
        return True, "Molecule is an N-acyl-alpha-amino acid (without specific stereochemistry)"

    return False, "Molecule is not an N-acyl-L-alpha-amino acid"