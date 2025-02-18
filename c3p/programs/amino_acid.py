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
        bool: True if molecule is classified as an amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid pattern (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O&H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for amino groups (-NH2, -NHR, -NR2)
    amino_group_pattern = Chem.MolFromSmarts("[NX3,NX4+][H,D1]")
    amino_group_matches = mol.GetSubstructMatches(amino_group_pattern)
    if len(amino_group_matches) < 1:
        return False, "No amino groups found"

    return True, "Contains carboxylic acid and amino group(s), classifiable as an amino acid"