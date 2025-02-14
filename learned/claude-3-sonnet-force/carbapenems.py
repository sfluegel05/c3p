"""
Classifies: CHEBI:46633 carbapenems
"""
"""
Classifies: CHEBI:49144 carbapenems
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    Carbapenems are a class of beta-lactam antibiotics with a carbapenem skeleton
    that is variously substituted and typically contains sulfur, amine, and/or hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbapenem, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carbapenem backbone pattern
    backbone_pattern = Chem.MolFromSmarts("[C@@H]1[C@H]2C(=O)N1C2")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No carbapenem backbone found"
    
    # Check for common functional groups
    sulfur_pattern = Chem.MolFromSmarts("S")
    amine_pattern = Chem.MolFromSmarts("N")
    hydroxy_pattern = Chem.MolFromSmarts("O[H]")
    
    has_sulfur = mol.HasSubstructMatch(sulfur_pattern)
    has_amine = mol.HasSubstructMatch(amine_pattern)
    has_hydroxy = mol.HasSubstructMatch(hydroxy_pattern)
    
    if not (has_sulfur or has_amine or has_hydroxy):
        return False, "Missing common functional groups (sulfur, amine, or hydroxyl)"
    
    # Additional checks for specific substituents or patterns
    # ...
    
    return True, "Contains carbapenem skeleton with common functional groups and substituents"