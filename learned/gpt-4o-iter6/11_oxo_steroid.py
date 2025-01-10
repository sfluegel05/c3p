"""
Classifies: CHEBI:47787 11-oxo steroid
"""
from rdkit import Chem

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an 11-oxo steroid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Attempting a more flexible pattern that captures the steroid backbone
    steroid_framework_pattern = Chem.MolFromSmarts("C1CC2CCCC3C=CC(=O)CCC3C2C1") 
    if steroid_framework_pattern is None or not mol.HasSubstructMatch(steroid_framework_pattern):
        return False, "No generalized steroid-like backbone found"
    
    # 11-oxo group: Refine to detect any typical position where a carbonyl (oxo) is adjacent in fused rings
    # Pattern allows flexibility in the actual carbon number where C=O is attached in a steroid-like framework
    oxo_11_generic_pattern = Chem.MolFromSmarts("C(=O)C1CCC2C(C1)CC3C4=C(CCC3C2)C(=O)C")
    if oxo_11_generic_pattern is None or not mol.HasSubstructMatch(oxo_11_generic_pattern):
        return False, "No 11-oxo group found in a suitable position for this framework"

    return True, "Contains steroid backbone with an oxo group likely at the 11th position"