"""
Classifies: CHEBI:33709 amino acid
"""
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

    # Look for carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Look for amino group (-NH2, -NH-, or -N-)
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(N=*)]")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"

    # If both patterns are found, it is an amino acid
    return True, "Contains at least one carboxylic acid group and one amino group"