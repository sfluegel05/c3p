"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
from rdkit import Chem

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    This class has a steroid backbone, a ketone at the third position,
    and alpha configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a relaxed steroid backbone pattern to capture structural diversity
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4CCCCC4C3CCC2C1") # This allows for more flexibility

    # Check for steroid backbone
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "Does not have a detectable steroid backbone"
    
    # Define the pattern for the 3-oxo group as a specific ketone
    three_oxo_pattern = Chem.MolFromSmarts("C(=O)C")  # Slightly relaxed to match variations

    # Check for the presence of a 3-oxo group
    if not mol.HasSubstructMatch(three_oxo_pattern):
        return False, "Missing the characteristic 3-oxo group"

    # Define specific pattern for 5alpha stereochemistry
    five_alpha_pattern = Chem.MolFromSmarts("[C@H]1CC[C@H]2[C@@H](C)C(=O)[C@@]3[C@H](C1)CC[C@@H]23")

    # Check for 5alpha configuration in the context of a steroid
    if not mol.HasSubstructMatch(five_alpha_pattern):
        return False, "Does not match 5alpha-steroid stereochemistry configuration"
    
    return True, "Matches 3-oxo-5alpha-steroid structure"