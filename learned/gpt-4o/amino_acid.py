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

    # Define the amino group pattern, capturing primary, secondary, and tertiary amines
    amino_patterns = [
        Chem.MolFromSmarts("[NX3;H2,H1]"),
        Chem.MolFromSmarts("[NX3;H0]")  # Tertiary or quaternary amines if necessary
    ]

    # Ensure the presence of an amino group
    if not any(mol.HasSubstructMatch(pattern) for pattern in amino_patterns):
        return False, "No amino group found"

    # Broaden patterns to capture potential amino acid-like backbones, including potential modifications or atypical arrangements
    broad_backbone_patterns = [
        Chem.MolFromSmarts("[NX3][C][C](=O)O"),
        Chem.MolFromSmarts("[NX3][C][CX4](=O)O"),
        Chem.MolFromSmarts("[NX3][C;H2][CX3](=O)O"),
        Chem.MolFromSmarts("[NX3][C;!R][CX3](=O)O")
    ]

    # Check for a broader sense of amino acid backbones to capture atypical modifications
    if not any(mol.HasSubstructMatch(pattern) for pattern in broad_backbone_patterns):
        return False, "No typical or modified amino acid backbone found"

    return True, "Contains carboxylic acid group and amino group, suitable for an amino acid structure"