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

    # Improved steroid backbone pattern: allows variability in ring positions
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CCC2[C@H]3[C@@H]4C5CC[C@H](C5)[C@@H]4[C@H]3CCC2C1")

    # Check for steroid backbone
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "Does not have a detectable steroid backbone"
    
    # Enhanced 3-oxo group detection, considering possible connected entities
    three_oxo_pattern = Chem.MolFromSmarts("C(=O)[C@]")

    # Check for the presence of a 3-oxo group
    if not mol.HasSubstructMatch(three_oxo_pattern):
        return False, "Missing the characteristic 3-oxo group"

    # Enhanced 5alpha stereochemistry pattern, ensuring position 5 is correctly matched
    five_alpha_stereochemistry = Chem.MolFromSmarts("[C@H](CC)[C@H](C(=O))")

    # Check for 5alpha configuration and that it's relevantly positioned
    if not mol.HasSubstructMatch(five_alpha_stereochemistry):
        return False, "Does not match 5alpha-steroid stereochemistry configuration"
    
    return True, "Matches 3-oxo-5alpha-steroid structure"