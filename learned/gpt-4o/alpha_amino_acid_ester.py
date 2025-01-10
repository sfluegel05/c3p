"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester typically has an amine (NH2) group on the alpha carbon of the carboxylic acid,
    which is further esterified.

    Args:
        smiles (str): SMILES string of the molecule
    Returns:
        bool: True if molecule is an alpha-amino acid ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for alpha-amino acid
    # Looking for an amine (N) attached to a carbon (C) which is adjacent to a carboxylic acid (COOH)
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[$([NX3][CX4][$([CX3](=O)[OX2H1])])]")

    # Define SMARTS pattern for ester linkage
    # Looking for an ester group O=C-O in the molecule
    ester_pattern = Chem.MolFromSmarts("[$([CX3](=O)[OX2])]")
    
    # Search for alpha-amino acid pattern
    if not mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return False, "No alpha-amino acid substructure found"
    
    # Search for ester group pattern
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"
    
    return True, "Contains an alpha-amino acid substructure with ester linkages"