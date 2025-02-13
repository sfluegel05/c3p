"""
Classifies: CHEBI:73754 thiosugar
"""
from rdkit import Chem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar is a carbohydrate derivative in which one or more of the oxygens
    or hydroxy groups of the parent sugar is replaced by sulfur or -SR groups.

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

    # General pattern for sugar backbone: flexible enough to cover various pyranoses and furanoses
    sugar_backbone_pattern = Chem.MolFromSmarts("[C@H]([C@H](O)*)[C@H](O)[C@H](O)[C@H](O)[C@@H]1O")  # More general template

    if not mol.HasSubstructMatch(sugar_backbone_pattern):
        return False, "No sugar backbone found"

    # Define sulfur substitution patterns
    sulfur_patterns = [
        Chem.MolFromSmarts("[#16]"),                   # Match sulfur
        Chem.MolFromSmarts("[C@H]S"),                  # S connected to a chiral center
        Chem.MolFromSmarts("O[SX4]"),                  # Sulfur replacing oxygen in pattern
        Chem.MolFromSmarts("C-S(=O)"),                 # Sulfur bonded to carbon, like sulfoxides
        Chem.MolFromSmarts("[SX2][C]"),                # Thioether configurations
    ]

    # Check for any sulfur substituting traditionally oxygen positions
    sulfur_substitution_found = any(mol.HasSubstructMatch(pat) for pat in sulfur_patterns)

    if not sulfur_substitution_found:
        return False, "No sulfur substitution found"
    
    return True, "Contains sugar backbone with sulfur substitution"