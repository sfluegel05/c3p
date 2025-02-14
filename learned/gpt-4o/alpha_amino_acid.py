"""
Classifies: CHEBI:33704 alpha-amino acid
"""
from rdkit import Chem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid has an amino group (NH2 or a derivative) on the carbon atom adjacent to (alpha to) the carboxylic acid (COOH or COO-) group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Alpha-amino acid generic patterns:
    # 1. [NX3][CX4](C)C(=O)O - matches typical alpha amino acids and considers charged states (-COO-)
    # 2. Include patterns for derivatives or more complex stereo-centers
    
    # Comprehensive alpha-amino pattern to cover common cases:
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[$([NX3][CX4H1R][CX3](=O)[O,H1,-1])]")

    if mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return True, "Contains alpha-amino acid structure"
    else:
        return False, "No alpha-amino acid pattern found"