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
    
    # Improved germacranolide patterns for variation in skeleton and functional groups

    # Pattern 1: Germacrane skeleton with a lactone and a 10-membered ring core
    germacranolide_pattern1 = Chem.MolFromSmarts("C1CCC2C(CCC2C1)C3=CC=CC(O3)=O")

    # Pattern 2: Germacrane core with ester linkage and macrocyclic lactone
    germacranolide_pattern2 = Chem.MolFromSmarts("C1C=C2CCC=C3C(C=CC3(C)O2)C1=O")

    # Pattern 3: Exocyclic double bonds, ester groups, and various substitutions
    germacranolide_pattern3 = Chem.MolFromSmarts("C1=CC=C2C=CC3C(C=CC3(C)O2)C1=O")
    
    # Check against multiple patterns
    if mol.HasSubstructMatch(germacranolide_pattern1):
        return True, "Contains germacranolide-like structure with lactone and 10-membered ring"
    if mol.HasSubstructMatch(germacranolide_pattern2):
        return True, "Contains germacranolide-like structure with ester linkage"
    if mol.HasSubstructMatch(germacranolide_pattern3):
        return True, "Contains germacranolide-like structure with exocyclic double bonds and substitutions"

    return False, "No identifiable germacranolide-like structure using current patterns"