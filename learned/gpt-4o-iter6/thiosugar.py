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
    
    # Detect the presence of a sugar moiety
    sugar_patterns = [
        Chem.MolFromSmarts("[C@@H]1O[C@H]([C@H](O)[C@@H](O)[C@H]1O)"), # pyranose
        Chem.MolFromSmarts("[C@@H]1O[C@H]([C@H](O)[C@H]1O)"),          # furanose
    ]
    
    has_sugar = any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns)
    if not has_sugar:
        return False, "No sugar backbone found"
    
    # Look for sulfur replacements (S or -SR)
    sulfur_patterns = [
        Chem.MolFromSmarts("[C-S]"),  # Any carbon-sulfur bond
        Chem.MolFromSmarts("[O-S]"),  # Oxygen replaced by sulfur
        Chem.MolFromSmarts("[S;D1]")  # Terminal sulfur
    ]
    
    has_sulfur = any(mol.HasSubstructMatch(pattern) for pattern in sulfur_patterns)
    if not has_sulfur:
        return False, "No sulfur replacement detected"
    
    return True, "Contains carbohydrate structure with sulfur substitution(s)"

# Example of how to use the function
# smiles_example = "CC(C)S[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O"  # isopropyl beta-D-thiogalactopyranoside
# print(is_thiosugar(smiles_example))