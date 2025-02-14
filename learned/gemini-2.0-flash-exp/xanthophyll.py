"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: CHEBI:27327 xanthophyll
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    Xanthophylls are oxygenated carotenoids with a polyene backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check number of carbons: should be at least 30
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 30:
        return False, "Too few carbon atoms for a xanthophyll."

    # Check for a longer polyene backbone, common for carotenoids
    # This is a simplified pattern, which matches C=C-C=C-C=C
    # It can be more complex for example "[CX3]=[CX3]-[CX3]([CH3])=[CX3]-[CX3]=[CX3]-[CX3]([CH3])=[CX3]"
    carotenoid_pattern = Chem.MolFromSmarts("[CX3]=[CX3]-[CX3]=[CX3]-[CX3]=[CX3]")
    if not mol.HasSubstructMatch(carotenoid_pattern):
      return False, "No basic carotenoid-like backbone structure found."

    # Check for total number of oxygen atoms: should be at least 1
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count == 0:
        return False, "No oxygen atoms present, not a xanthophyll."
    
    # Check for multiple oxygen-containing functional groups (hydroxyl, carbonyl, epoxide)
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    carbonyl_pattern = Chem.MolFromSmarts("C=[OX1]")
    epoxide_pattern = Chem.MolFromSmarts("C1OC1")
    oxygenated_groups = mol.GetSubstructMatches(hydroxy_pattern) + mol.GetSubstructMatches(carbonyl_pattern) + mol.GetSubstructMatches(epoxide_pattern)
    if len(oxygenated_groups) < 1:
        return False, "No characteristic oxygen-containing functional groups found."
    if oxygen_count < 2:
      return False, "Xanthophylls usually contain multiple oxygen atoms."

    # Check for glycoside like structures using a slightly more permissive SMARTS
    glycoside_pattern = Chem.MolFromSmarts("OC[C@H]1[C@@H](O)[C@@H](O)[C@H]([OX2])[C@@H](CO)O[C@H]1")
    if mol.HasSubstructMatch(glycoside_pattern):
      return True, "Contains a carotenoid-like backbone with characteristic oxygen-containing functional groups including glycosides."
      
    return True, "Contains a carotenoid-like backbone with oxygen-containing functional groups."