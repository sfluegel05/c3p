"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester is formed by the formal condensation of an alpha-amino acid with an alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid ester, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into a molecular structure
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for an alpha-amino acid ester
    # [C](=O)O[C] for ester linkage, with an amino group `[NH2,NH]` on the alpha-carbon
    alpha_amino_ester_pattern = "[C;H1,H2,H3]([CH1,CH2]([NH2,NH])[C](=O)O[C])"

    # Convert the pattern to a molecule
    pattern_mol = Chem.MolFromSmarts(alpha_amino_ester_pattern)
    if pattern_mol is None:
        return None, None  # There was an error creating the pattern

    # Check if the molecule matches the alpha-amino acid ester pattern
    if mol.HasSubstructMatch(pattern_mol):
        return True, "Contains an alpha-amino acid backbone with ester linkage"

    return False, "Does not match the alpha-amino acid ester pattern"