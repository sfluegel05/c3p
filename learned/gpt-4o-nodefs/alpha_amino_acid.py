"""
Classifies: CHEBI:33704 alpha-amino acid
"""
from rdkit import Chem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid typically has a central alpha carbon atom bonded to an amino group (NH2),
    a carboxylic acid group (COOH or COO-), and a side chain.

    Args:
        smiles (str): SMILES string of the molecule
    Returns:
        bool: True if molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for alpha amino acids: a central carbon bonded to N, C(=O)[O,H-], and side chains (excluding hydrogens)
    # This pattern is flexible to capture common amino acids as well as some modified residues
    alpha_amino_pattern = Chem.MolFromSmarts("[C;H1,H2,H3](N)([CX3](=O)[OX1H0-,OX2H1])")

    if mol.HasSubstructMatch(alpha_amino_pattern):
        return True, "Structure matches an alpha-amino acid pattern"

    return False, "No match to the alpha-amino acid structure"