"""
Classifies: CHEBI:73011 germacranolide
"""
from rdkit import Chem

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide based on its SMILES string.
    A germacranolide is a sesquiterpene lactone based on a germacrane skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a germacranolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Refined SMARTS pattern for a germacranolide (fused bicyclic + lactone)
    # This pattern captures:
    # - A bicyclic decalin-like core: Two fused six-membered rings
    # - A five-membered lactone ring fused or attached to the core
    germacranolide_pattern = Chem.MolFromSmarts("C1=CCC2C(C=CC2O1)CC(=O)O")
    
    # Check for germacranolide structure
    if not mol.HasSubstructMatch(germacranolide_pattern):
        return False, "No germacranolide-like structure found"
    
    return True, "Contains a bicyclic germacrane skeleton with an embedded lactone group"