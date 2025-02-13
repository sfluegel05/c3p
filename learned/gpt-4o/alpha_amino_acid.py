"""
Classifies: CHEBI:33704 alpha-amino acid
"""
from rdkit import Chem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid has an amino group on the carbon atom adjacent to (alpha to) the carboxylic acid group.

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
    
    # Improved pattern recognizing diverse alpha-amino acid structures:
    # Account for stereochemistry, ring structures, and complex side chains
    patterns = [
        Chem.MolFromSmarts("[$([NX3R][C@@H1R][CX3](=O)[OX1-,OX2H1])]"),  # Stereochemistry + typical alpha-amino (COOH/COO-)
        Chem.MolFromSmarts("[CX3H1](N)[CX3](=O)[OX1-,OX2H1]"),            # Simple backbone with variations
        Chem.MolFromSmarts("[NX3H,NX4H2]-[CX4H]([*])[CX3](=O)[OX1-,OX2H1]"),  # Generic flexibility for various configurations
    ]
    
    # Check for matches with any of the defined patterns
    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains alpha-amino acid structure"

    return False, "No alpha-amino acid pattern found"