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
    Carbapenems are a class of beta-lactam antibiotics that have a carbapenem skeleton,
    which consists of a four-membered beta-lactam ring fused to a five-membered ring
    with a sulfur atom and a carboxylate group attached.

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
    backbone_pattern = Chem.MolFromSmarts("[C@@H]1[C@H]2C(=O)N1C2SC3=C(N4C(=O)[C@@]3([C@H]([C@@H]4C)O)C)C(O)=O")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No carbapenem backbone found"
    
    # Check for common functional groups and substituents
    hydroxy_pattern = Chem.MolFromSmarts("O[H]")
    amine_pattern = Chem.MolFromSmarts("N")
    
    has_hydroxy = mol.HasSubstructMatch(hydroxy_pattern)
    has_amine = mol.HasSubstructMatch(amine_pattern)
    
    if not (has_hydroxy or has_amine):
        return False, "Missing common functional groups (hydroxyl or amine)"
    
    return True, "Contains carbapenem skeleton with common functional groups and substituents"