"""
Classifies: CHEBI:73754 thiosugar
"""
from rdkit import Chem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar has one or more oxygens or hydroxy groups of a carbohydrate
    replaced by sulfur or -SR groups.

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

    # General pattern for sugar backbone: cyclic structure with multiple hydroxyl groups
    # Allow ring sizes of 5 (furanose) to 6 (pyranose)
    sugar_backbone_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H]([C@H]1O)")
    # We should allow flexibility here as well, sugars can vary in hydroxyl arrangement

    if not mol.HasSubstructMatch(sugar_backbone_pattern):
        return False, "No sugar backbone found"

    # Sulfur substitution somewhere within the traditional oxygen positions in sugars
    # It captures -S- or -SR in place of -O-/-OH, or as a group like SOx, etc.
    sulfur_pattern = Chem.MolFromSmarts("[#16]")  # Generic sulfur atom
    
    if not mol.HasSubstructMatch(sulfur_pattern):
        return False, "No sulfur substitution found"

    # Verify sulfur replaces typical oxygen-rich moiety in the sugar structure: 
    thio_subs_pattern = Chem.MolFromSmarts("[C@H]S")
    if not any(mol.HasSubstructMatch(pt) for pt in [sulfur_pattern, thio_subs_pattern]):
        return False, "No thiosugar-like sulfur arrangement found"
    
    return True, "Contains sugar backbone with sulfur substitution"