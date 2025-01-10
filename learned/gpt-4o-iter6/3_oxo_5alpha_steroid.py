"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
from rdkit import Chem

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    This class of compounds has a steroid backbone, a ketone at the third position, and alpha configuration at position 5.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define 3-oxo group pattern (ketone at position 3, generic for alpha and beta stereochemistry)
    three_oxo_pattern = Chem.MolFromSmarts("C(=O)[C@@H]")  # Adjust (add more valid patterns)
    
    # Check for 3-oxo group
    if not mol.HasSubstructMatch(three_oxo_pattern):
        return False, "Missing 3-oxo (ketone) group"
    
    # Steroid backbone pattern with 5alpha configuration
    five_alpha_steroid_pattern = Chem.MolFromSmarts("[C@H]1CC[C@@H]2[C@H]3CC[C@H]4[C@@H](C(=O))CC[C@]4(C)[C@@H]3CC[C@]12C")
    
    # Check for 5alpha configuration and steroid-like structure
    if not mol.HasSubstructMatch(five_alpha_steroid_pattern):
        return False, "Does not match 5alpha-steroid configuration"
    
    return True, "Matches 3-oxo-5alpha-steroid structure"