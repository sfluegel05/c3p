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
    
    # Parse the SMILES string into a molecular object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the carboxylic acid pattern including the ionized and non-ionized forms
    carboxylic_patterns = [
        Chem.MolFromSmarts("C(=O)[O-]"),  # Ionized form
        Chem.MolFromSmarts("C(=O)O")      # Non-ionized form
    ]
    
    # Check if any carboxylic pattern exists
    if not any(mol.HasSubstructMatch(pattern) for pattern in carboxylic_patterns):
        return False, "No carboxylic acid group found"

    # Define the amino group pattern, capturing primary and secondary amines
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1]")

    # Ensure the presence of an amino group
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"

    # Define patterns that capture typical amino acid backbone connectivity
    # N-C-C(-C)-C(=O)O / N-C-C(-C)-C(=O)[O-]
    backbone_patterns = [
        Chem.MolFromSmarts("[NX3][CH2][CH](C)[C](=O)O"),
        Chem.MolFromSmarts("[NX3][CH2][CH](C)[C](=O)[O-]"),
        Chem.MolFromSmarts("[NX3][CH2][CH2][C](=O)[O-]")
    ]
    
    # Check that at least one backbone pattern is present
    if not any(mol.HasSubstructMatch(pattern) for pattern in backbone_patterns):
        return False, "No typical amino acid backbone found"

    return True, "Contains carboxylic acid group and amino group in a typical amino acid structure"