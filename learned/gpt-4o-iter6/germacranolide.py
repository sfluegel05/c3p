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
    
    # Expanded germacranolide patterns for variation in skeleton and functional groups

    # Pattern 1: Basic germacrane skeleton with lactone
    germacranolide_pattern1 = Chem.MolFromSmarts("C1CC2CCC(C1)C(C)CC2OC(=O)")

    # Pattern 2: Alternate skeleton configurations (larger macrocyclic structures)
    germacranolide_pattern2 = Chem.MolFromSmarts("O=C1OC2C=CC=CCC[C@@H]2CCC1")

    # Pattern 3: Potential for ester substitutions and exocyclic double bonds
    germacranolide_pattern3 = Chem.MolFromSmarts("C1[C@H]2CC(C1)C=CC[C@H](C2)C(=O)O")
    
    # Check against multiple patterns
    if mol.HasSubstructMatch(germacranolide_pattern1):
        return True, "Contains germacranolide-like structure with basic germacrane skeleton and lactone"
    if mol.HasSubstructMatch(germacranolide_pattern2):
        return True, "Contains germacranolide-like structure with alternate germacranolide skeleton"
    if mol.HasSubstructMatch(germacranolide_pattern3):
        return True, "Contains germacranolide-like structure with ester substitutions"

    return False, "No identifiable germacranolide-like structure using current patterns"