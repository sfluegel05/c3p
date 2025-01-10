"""
Classifies: CHEBI:73754 thiosugar
"""
from rdkit import Chem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar is a carbohydrate derivative in which one or more of the oxygens
    or hydroxy groups of the parent carbohydrate is replaced by sulfur or -SR.

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
    
    # Detect the presence of a sugar moiety using an expanded set of patterns
    sugar_patterns = [
        Chem.MolFromSmarts("[C@@H]1O[C@H]([C@@H](O)*)[C@@H](O)C[C@H]1O"),  # pyranose variant
        Chem.MolFromSmarts("[C@@H]1O[C@H]1"),                             # general sugar cycle
        Chem.MolFromSmarts("[O-C@[C@](O)(C-O)C-O]")                       # simplified open chain
    ]
    
    has_sugar = any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns)
    if not has_sugar:
        return False, "No sugar backbone found"
    
    # Look for sulfur replacements in positions where sugars usually have oxygens
    sulfur_patterns = [
        Chem.MolFromSmarts("[O;!H1]-S"),  # oxygens directly replaced by sulfur
        Chem.MolFromSmarts("[C-S]"),      # carbon-sulfur bonds adjacent to oxy-like positions
        Chem.MolFromSmarts("[S;X2]")      # sulfur in right oxidation state for substitution
    ]
    
    valid_sulfur_replacement = any(mol.HasSubstructMatch(pattern) for pattern in sulfur_patterns)
    if not valid_sulfur_replacement:
        return False, "No valid sulfur replacement detected"
    
    return True, "Contains sugar structure with sulfur substitution(s)"