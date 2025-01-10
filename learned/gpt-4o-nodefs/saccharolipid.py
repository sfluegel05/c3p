"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    It checks for features typical of saccharolipids: fatty acid chains, sugar moieties, and ester/ether linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Saccharolipids typically include long hydrocarbon chains (fatty acids)
    long_chain_pattern = Chem.MolFromSmarts("C(CCCC)(CCCC)(CCCC)")  # Simplified example for a long alkyl chain
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long hydrocarbon chains (saccharolipid characteristic) detected"
    
    # Look for sugar moieties - attempt to capture common sugar patterns broadly
    sugar_patterns = [
        Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](CO)O1"),  # Glucose-like
        Chem.MolFromSmarts("O[C@H]1CO[C@@H](O)[C@H](O)[C@H](O)[C@H]1O")]  # Another cyclic sugar
    if not any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns):
        return False, "No sugar moieties detected"
    
    # Check for ester or glycosidic linkages
    linkage_patterns = [
        Chem.MolFromSmarts("[C](=O)[O][C]"),  # Ester
        Chem.MolFromSmarts("O[C@H]")]  # Glycosidic
    if not any(mol.HasSubstructMatch(pattern) for pattern in linkage_patterns):
        return False, "No ester or glycosidic linkages found"

    return True, "Contains features typical of saccharolipids: fatty acids, sugars, and characteristic linkages"