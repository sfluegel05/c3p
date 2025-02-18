"""
Classifies: CHEBI:33709 amino acid
"""
from rdkit import Chem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid is defined as having a carboxylic acid group and an amino group,
    both directly connected or through non-hindering atoms/bonds, to the alpha carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classifiable as an amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for amino acid backbone
    # [N] represents the nitrogen (part of the amino group), 
    # [C](=O)O represents the carboxylic acid
    # X represents any atom that might be between the nitrogen and the alpha carbon.
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4H]([CX3](=O)[OX2H1])")

    if mol.HasSubstructMatch(amino_acid_pattern):
        return True, "Contains amino acid backbone structure"

    return False, "No amino acid backbone structure found"