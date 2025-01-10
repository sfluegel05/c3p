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
        bool: True if the molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the steroid backbone pattern generally as a tetracyclic steroid structure
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CCC2C3C4CCCCC4C3CCC2C1")  # More typical steroid backbone

    # Check for steroid backbone
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "Does not have a steroid backbone"
    
    # Define 3-oxo group pattern: a ketone group consistent with steroid structure
    three_oxo_pattern = Chem.MolFromSmarts("C[C](=O)[C@H]")  
    
    # Check for 3-oxo group
    if not mol.HasSubstructMatch(three_oxo_pattern):
        return False, "Missing 3-oxo (ketone) group at expected position"
    
    # Define specific pattern for 5alpha stereochemistry
    # This pattern assumes that the C-5 should be alpha to the plane and appropriately stereocenters are identified
    five_alpha_pattern = Chem.MolFromSmarts("[C@@H]1CC[C@H]2[C@@H](C=C3[C@H]CC[C@@H]23)C1")

    # Check for 5alpha configuration in the context of a steroid
    if not mol.HasSubstructMatch(five_alpha_pattern):
        return False, "Does not match 5alpha-steroid stereochemistry configuration"
    
    return True, "Matches 3-oxo-5alpha-steroid structure"