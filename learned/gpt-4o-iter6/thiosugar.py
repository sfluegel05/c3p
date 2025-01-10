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
    
    # Detect the presence of a sugar moiety (cycles with oxygen atoms and carbon pattern)
    sugar_patterns = [
        Chem.MolFromSmarts("[C@@H]1O[C@H]([C@@H](O)[C@@H](O)[C@@H]1O)"),  # pyranose 6-membered ring
        Chem.MolFromSmarts("[C@@H]1O[C@H]([C@@H](O)[C@@H]1O)"),           # furanose 5-membered ring
        Chem.MolFromSmarts("[C-O-C-O-C-O]")                              # general sugar pattern
    ]
    
    has_sugar = any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns)
    if not has_sugar:
        return False, "No sugar backbone found"
    
    # Look for sulfur replacements (S or -SR) where sulfur is bound to a carbon adjacent to the sugar
    sulfur_patterns = [
        Chem.MolFromSmarts("[C-S]"),  # Any carbon-sulfur bond
        Chem.MolFromSmarts("[O-S]"),  # Direct replacement of oxygen
        Chem.MolFromSmarts("[S;D1]")  # Terminal sulfur 
    ]
    
    valid_sulfur_replacement = False
    for pattern in sulfur_patterns:
        if mol.HasSubstructMatch(pattern):
            valid_sulfur_replacement = True
            break

    if not valid_sulfur_replacement:
        return False, "No valid sulfur replacement detected"
    
    return True, "Contains sugar structure with sulfur substitution(s)"