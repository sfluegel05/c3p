"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid contains both a carbohydrate moiety and a lipid component.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define broader patterns for lipid and carbohydrate components
    # Lipid pattern: optional functional groups around long carbon chains
    lipid_pattern = Chem.MolFromSmarts("C(=O)OC([C@@H,CH,C])([C@@H,CH,C])C")  # Esterified long chain
    # Carbohydrate pattern: generic sugar rings, considering anomeric centers
    carb_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1")
    
    # Check for presence of lipid component
    if not mol.HasSubstructMatch(lipid_pattern):
        return False, "No adequate long hydrocarbon chains (lipid component) found"
    
    # Check for presence of carbohydrate component
    if not mol.HasSubstructMatch(carb_pattern):
        return False, "No adequate carbohydrate moiety found"
    
    return True, "Contains both lipid and carbohydrate components indicating a saccharolipid"