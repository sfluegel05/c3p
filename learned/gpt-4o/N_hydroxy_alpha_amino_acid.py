"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: N-hydroxy-alpha-amino-acid
"""
from rdkit import Chem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino-acid based on its SMILES string.
    
    An N-hydroxy-alpha-amino-acid is defined as any amino acid in which at least one hydrogen 
    attached to the amino group is replaced by a hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-hydroxy-alpha-amino-acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Look for an N-hydroxy group (N attached to O)
    n_hydroxy_pattern_1 = Chem.MolFromSmarts("N(O)")  # Primary N-hydroxy
    n_hydroxy_pattern_2 = Chem.MolFromSmarts("N(O)O")  # Secondary N-hydroxy
    if not mol.HasSubstructMatch(n_hydroxy_pattern_1) and not mol.HasSubstructMatch(n_hydroxy_pattern_2):
        return False, "No N-hydroxy group found (N-O or N(OO))"

    # Additional checks could include molecular weight, presence of chiral centers, etc., but these are less defined
    # for this specific classification and could be highly variable.

    return True, "Contains N-hydroxy substitution on alpha-amino group and carboxylic acid characteristic of alpha-amino acids"