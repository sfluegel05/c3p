"""
Classifies: CHEBI:28963 amino sugar
"""
from rdkit import Chem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is a sugar having one or more alcoholic hydroxy groups replaced
    by substituted or unsubstituted amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flexible pattern for hexose sugars
    # This pattern should capture a six-membered ring with hydroxyl groups
    sugar_pattern = Chem.MolFromSmarts("C1(CO)OC(O)C(O)C(O)C1")

    # Flexible pattern for furanose sugars
    sugar_pattern_furanose = Chem.MolFromSmarts("C1(CO)OC(O)C(O)C1")

    # Check if there's a sugar backbone
    if mol.HasSubstructMatch(sugar_pattern) or mol.HasSubstructMatch(sugar_pattern_furanose):
        # Define a pattern for amino substitution
        amino_pattern = Chem.MolFromSmarts("[NX3,NX4][C;!$(C=O)]")  # amino group attached to carbon not a part of carbonyl
        if mol.HasSubstructMatch(amino_pattern):
            return True, "Contains sugar backbone with hydroxyl groups replaced by amino groups"
        else:
            return False, "Sugar backbone found but no amino substitution detected"
    
    return False, "No sugar backbone found"