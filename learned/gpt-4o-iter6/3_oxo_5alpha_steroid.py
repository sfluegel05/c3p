"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
from rdkit import Chem

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    This class of compounds has a steroid backbone, a ketone at the third position, 
    and alpha configuration at position 5.
    
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
    
    # Define 3-oxo group pattern: General ketone group
    three_oxo_pattern = Chem.MolFromSmarts("C(=O)C")  
    
    # Check for 3-oxo group
    if not mol.HasSubstructMatch(three_oxo_pattern):
        return False, "Missing 3-oxo (ketone) group"
    
    # Steroid backbone pattern allowing some flexibility in the rings
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4C(CCC4C3)C2C1")  # General steroid scaffold
    
    # Check for steroid backbone
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "Does not match steroid-like structure"
    
    # Specific pattern for 5alpha stereochemistry (allow more stereochemical flexibility)
    five_alpha_pattern = Chem.MolFromSmarts("[C@@H]([C@@H]1)[C@@H]2[C@H]3CC[C@@H](C2)C[C@H]1")
    
    # Check for 5alpha configuration in the context of a steroid
    if not mol.HasSubstructMatch(five_alpha_pattern):
        return False, "Does not match 5alpha-steroid configuration"
    
    return True, "Matches 3-oxo-5alpha-steroid structure"