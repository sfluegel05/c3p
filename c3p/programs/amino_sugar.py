"""
Classifies: CHEBI:28963 amino sugar
"""
from rdkit import Chem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is any sugar having one or more alcoholic hydroxy groups
    replaced by substituted or unsubstituted amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for standard C-O-C or C-C-O patterns found in sugars
    sugar_pattern = Chem.MolFromSmarts("C(O)C(O)C(O)")  # Simplified sugar-like structure
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar-like backbone found"

    # Look for NH2 or substituted amino group replacing an OH
    amino_pattern = Chem.MolFromSmarts("[CX4][NH2,NH,R,NR2,C(=O)NR]")  # Simplified amino pattern
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group replacing hydroxy group found"

    return True, "Contains sugar-like backbone with one or more hydroxyl groups replaced by amino groups"