"""
Classifies: CHEBI:23437 cyanohydrin
"""
"""
Classifies: CHEBI:39142 cyanohydrin
"""
from rdkit import Chem

def is_cyanohydrin(smiles: str):
    """
    Determines if a molecule is a cyanohydrin based on its SMILES string.
    A cyanohydrin is an alpha-hydroxynitrile resulting from the addition of HCN to an aldehyde or ketone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyanohydrin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pattern 1: Hydroxyl and nitrile on the same carbon
    pattern1 = Chem.MolFromSmarts("[C]([OH])[C]#[N]")
    # Pattern 2: Hydroxyl on adjacent carbon to nitrile
    pattern2 = Chem.MolFromSmarts("[C][C]([OH])#[N]")
    
    # Check for either pattern
    has_pattern1 = mol.HasSubstructMatch(pattern1)
    has_pattern2 = mol.HasSubstructMatch(pattern2)
    
    if has_pattern1 or has_pattern2:
        return True, "Contains hydroxyl and nitrile groups on adjacent carbons or the same carbon"
    
    return False, "No hydroxyl and nitrile groups in required positions"