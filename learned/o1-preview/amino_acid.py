"""
Classifies: CHEBI:33709 amino acid
"""
"""
Classifies: amino acid
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

    # Carboxylic acid group: matches both protonated and deprotonated forms
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Amino group: primary, secondary, or protonated amines (excluding amides and imines)
    amino_group_pattern = Chem.MolFromSmarts("[N;!$(N-C=O);!$(N-[*]=[*]);H1,H2,H3;+0]")
    if not mol.HasSubstructMatch(amino_group_pattern):
        return False, "No amino group found"

    return True, "Contains carboxylic acid group and one or more amino groups"