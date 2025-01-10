"""
Classifies: CHEBI:17002 cholesteryl ester
"""
from rdkit import Chem

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    A cholesteryl ester is defined by esterification of the 3-hydroxy group of cholesterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cholesteryl ester, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the general sterol core pattern including cholesteryl specifics (as flexible as needed)
    cholesterol_backbone_pattern = Chem.MolFromSmarts("C[C@H]1CC[C@]2(C)[C@H]3CC[C@H]4[C@H](CC[C@]4(C)[C@@H]3CC=C2C1)O")
    if not mol.HasSubstructMatch(cholesterol_backbone_pattern):
        return False, "No cholesterol steroid backbone found"
    
    # Define pattern for ester linkage
    ester_linkage_pattern = Chem.MolFromSmarts("C(=O)O[C@H]1")
    if not mol.HasSubstructMatch(ester_linkage_pattern):
        return False, "No ester linkage found"

    # Find if the ester linkage is correctly connected to the cholesterol structure
    ester_matches = mol.GetSubstructMatches(ester_linkage_pattern)
    cholesterol_matches = mol.GetSubstructMatches(cholesterol_backbone_pattern)

    # Cross-check the specific connection point for a more precise ester linkage connection
    if not any(set(ester_match) & set(cholesterol_match) for ester_match in ester_matches for cholesterol_match in cholesterol_matches):
        return False, "Ester linkage not specifically with cholesterol's 3-hydroxy position"

    return True, "Contains cholesterol backbone with specific ester linkage indicative of cholesteryl ester"