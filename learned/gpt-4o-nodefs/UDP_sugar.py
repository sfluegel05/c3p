"""
Classifies: CHEBI:17297 UDP-sugar
"""
from rdkit import Chem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a UDP-sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Modified SMARTS patterns for components of UDP-sugar
    # General Uracil-like moiety (simplifying to account for potential tautomeric forms)
    uracil_pattern = Chem.MolFromSmarts("c1cc(=O)[nH]c(=O)n1") 
    
    # Ribose ring with possible diphosphate attachment
    diphosphate_ribose_pattern = Chem.MolFromSmarts("O[P](=O)(O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O") 
    
    # Look for uracil presence (allow flexibility for tautomers and variations in representation)
    if not mol.HasSubstructMatch(uracil_pattern):
        return False, "No uracil moiety found"

    # Check for diphosphate linkage which may include or exclude a specific ribose configuration
    if not mol.HasSubstructMatch(diphosphate_ribose_pattern):
        return False, "No compatible diphosphate ribose linkage found"
    
    # General sugar moiety linked via glycosidic bond, now including variability for sugar loop size
    sugar_glycosidic_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@H](O*)[C@@H]1O") 
    if not mol.HasSubstructMatch(sugar_glycosidic_pattern):
        return False, "No suitable sugar moiety linked via glycosidic bond"

    return True, "Contains uracil, compatible diphosphate ribose linkage, and appropriate sugar moiety characteristic of UDP-sugars"