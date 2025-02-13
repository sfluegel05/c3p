"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: CHEBI:36675 saccharolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is a lipid that contains a carbohydrate moiety.

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
    
    # Check for the presence of a carbohydrate substructure
    carbohydrate_pattern = Chem.MolFromSmarts("OC[C@H]([C@H](O)[C@H](O)[C@H](O)O)O")
    if not mol.HasSubstructMatch(carbohydrate_pattern):
        return False, "No carbohydrate moiety found"
    
    # Check for the presence of a lipid substructure
    lipid_pattern = Chem.MolFromSmarts("[C;H3][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2]")
    if not mol.HasSubstructMatch(lipid_pattern):
        return False, "No lipid moiety found"
    
    # Count the number of carbohydrate and lipid moieties
    carbohydrate_count = len(mol.GetSubstructMatches(carbohydrate_pattern))
    lipid_count = len(mol.GetSubstructMatches(lipid_pattern))
    
    if carbohydrate_count > 0 and lipid_count > 0:
        return True, f"Contains {carbohydrate_count} carbohydrate moiety(ies) and {lipid_count} lipid moiety(ies)"
    else:
        return False, "Does not meet the criteria for a saccharolipid"