"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
from rdkit import Chem

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    3-oxo-5alpha-steroids have a steroid backbone, a ketone at the third position, and alpha configuration at position 5.
    
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
    
    # Define 3-oxo pattern (ketone group at position 3)
    # This is a simplified pattern, assuming numeric position information or more context may be needed in a real scenario
    three_oxo_pattern = Chem.MolFromSmarts("C(=O)[C@H]")  # C=O and adjacent to a chiral center
    
    # Check for 3-oxo group
    if not mol.HasSubstructMatch(three_oxo_pattern):
        return False, "Missing 3-oxo (ketone) group"
    
    # Define 5alpha-steroid pattern (5alpha configuration)
    # Assuming a generic steroid backbone structure plus specific chiral centers
    five_alpha_pattern = Chem.MolFromSmarts("[C@@H]1(CC[C@@H]([C@H]2CC[C@@H]3[C@H](CC[C@H]23)C1)C(=O))")
    
    # Check for 5alpha configuration within a steroid-like structure
    if not mol.HasSubstructMatch(five_alpha_pattern):
        return False, "Does not match 5alpha-steroid configuration"
    
    # If more detailed checking of steroid backbone was needed, that could be added here
    
    return True, "Matches 3-oxo-5alpha-steroid structure"