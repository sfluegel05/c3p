"""
Classifies: CHEBI:33709 amino acid
"""
from rdkit import Chem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid is defined as a carboxylic acid containing one or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (C(=O)O)
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)[O-]")  # Slightly modified for flexibility
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"

    # Check for amino group (N)
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0]")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"

    # Check for typical amino acid structure (free or side-chain modifications)
    amino_acid_patterns = [
        Chem.MolFromSmarts("[C;X4](N)([C;X4])(C(=O)[O-])"),  # General alpha amino acids
        Chem.MolFromSmarts("[C;X4](C(=O)[O-])([N;X3])"),     # Beta or gamma alternative structure
        Chem.MolFromSmarts("[NX3;H2,H1,H0]-[C;X4](C=O)[O-]") # Variations with N-C-C-C backbone
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in amino_acid_patterns):
        return False, "No suitable amino acid connectivity pattern found"

    return True, "Contains carboxylic acid group and amino group in a typical amino acid structure"