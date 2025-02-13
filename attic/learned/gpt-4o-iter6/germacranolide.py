"""
Classifies: CHEBI:73011 germacranolide
"""
from rdkit import Chem

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    A germacranolide is a sesquiterpene lactone based on the germacrane skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a germacranolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for germacrane and lactone structures
    germacrane_pattern = Chem.MolFromSmarts("C1=CC=C(C=C1)C2CCCCC2")  # Simplified pattern for a decalin structure
    lactone_pattern = Chem.MolFromSmarts("C1OC(=O)C=C1")  # Simple lactone 5-membered ring
    
    # Check for germacrane-like structure
    if not mol.HasSubstructMatch(germacrane_pattern):
        return False, "No germacrane-like structure found"
        
    # Check for lactone group
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone group found"
    
    return True, "Contains a germacrane-like skeleton with an embedded lactone group"