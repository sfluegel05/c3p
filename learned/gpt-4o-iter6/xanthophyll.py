"""
Classifies: CHEBI:27325 xanthophyll
"""
from rdkit import Chem

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    A xanthophyll is a carotenoid that is oxygenated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # General pattern of polyene chain characteristic in carotenoids
    polyene_pattern = Chem.MolFromSmarts("[#6]=[#6](-[#6]=[#6])*,=,=[#6]")
    if polyene_pattern and not mol.HasSubstructMatch(polyene_pattern):
        return False, "Basic polyene chain character of carotenoids not found"
    
    # Check for the presence of any oxygen atom
    oxygen_presence = any(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())
    if not oxygen_presence:
        return False, "No oxygen detected; required for xanthophyll"

    # Check for potentially missed oxygen functionalities
    possible_oxy_func_groups = ["[OX2H]", "[CX3]=[OX1]", "[CX3][OX2R]"]
    for pattern in possible_oxy_func_groups:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return True, "Contains features of an oxygenated carotenoid (xanthophyll)"
    
    return False, "Lacks sufficient features to be classified as a xanthophyll"