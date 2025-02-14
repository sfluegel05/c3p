"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
from rdkit import Chem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for alpha-amino acid pattern
    # Carbon with an attached NH2 and COOH, and verify L-configuration
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[C@@H](N)C(=O)O")
    if not mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return False, "Does not match the L-alpha-amino acid structure"

    # The molecule has the required pattern for an L-alpha-amino acid
    return True, "Contains L-alpha-amino acid structure"