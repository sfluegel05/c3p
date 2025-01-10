"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: thiosugar
"""
from rdkit import Chem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar has a carbohydrate backbone with one or more oxygens or hydroxy
    groups replaced by sulfur or -SR groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for a common carbohydrate ring pattern
    # Simplified match for sugar-like patterns (5- or 6-membered rings with hydroxyl groups)
    carbohydrate_ring_pattern_5 = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@@H]1")
    carbohydrate_ring_pattern_6 = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@@H]1")
    
    if not (mol.HasSubstructMatch(carbohydrate_ring_pattern_5) or mol.HasSubstructMatch(carbohydrate_ring_pattern_6)):
        return False, "No carbohydrate-like ring structure found"
    
    # Check for sulfur substitutions
    sulfur_pattern = Chem.MolFromSmarts("[SH,X2]")
    sulfur_match = mol.HasSubstructMatch(sulfur_pattern)
    
    if sulfur_match:
        return True, "Sulfur found as substitution in carbohydrate-like structure"

    return False, "No sulfur substitution found in expected carbohydrate structure locations"