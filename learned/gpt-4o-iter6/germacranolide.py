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
    
    # Refined pattern: a decalin-like motif, a lactone ring, and flexible connections
    # We can try to match:
    # - A flexible terpene-like skeleton: Possible variants of bicyclic and lactone-containing motifs
    # For simplicity, an initial SMARTS including a germacrane with potential lactones:
    germacranolide_pattern1 = Chem.MolFromSmarts("C1CC2CCC(C1)C(C)CC2OC(=O)[C;R2]")  # A basic example germacrane-lactone
    
    # You might need multiple patterns due to variations in structure
    # Consider also additional lactone-containing patterns
    
    # Check for germacrane-like skeleton and lactone group
    if mol.HasSubstructMatch(germacranolide_pattern1):
        return True, "Contains germacranolide-like structure with terpene skeleton and lactone group"
    
    return False, "No identifiable germacranolide-like structure using current patterns"