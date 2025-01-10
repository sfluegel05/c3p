"""
Classifies: CHEBI:73754 thiosugar
"""
from rdkit import Chem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar is a carbohydrate derivative in which one or more oxygens
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
    
    # Detect the presence of a sugar moiety with a specific C-O pattern for known sugars
    try:
        sugar_patterns = [
            Chem.MolFromSmarts("[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@@H](O1)"),  # pyranosyl
            Chem.MolFromSmarts("[C@H]1O[C@@H](O)[C@H](O)[C@H]1"),              # furanosyl
            Chem.MolFromSmarts("C1OC(CO)C(O)C1"),                             # ring-like
        ]
    except:
        return None, None
    
    has_sugar = any(pattern is not None and mol.HasSubstructMatch(pattern) for pattern in sugar_patterns)
    if not has_sugar:
        return False, "No sugar backbone found"
    
    # Look for sulfur replacements in positions where sugars usually have oxygens
    try:
        sulfur_patterns = [
            Chem.MolFromSmarts("[O,S][C;X4]1[O,S][C;X4]([O,S])[C;X4]([O,S])[C;X4]1[O,S]"),  # pyranosyl with S
            Chem.MolFromSmarts("[C@H](O)[C@H](S)[C@H](O)"),                                 # S substitution
        ]
    except:
        return None, None
    
    valid_sulfur_replacement = any(pattern is not None and mol.HasSubstructMatch(pattern) for pattern in sulfur_patterns)
    if not valid_sulfur_replacement:
        return False, "No valid sulfur replacement detected"

    return True, "Contains sugar structure with sulfur substitution(s)"