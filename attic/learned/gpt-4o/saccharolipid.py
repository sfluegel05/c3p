"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is characterized by the presence of both a lipid and a carbohydrate moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define broader patterns for carbohydrate moieties (sugar rings or fragments)
    carbohydrate_patterns = [
        Chem.MolFromSmarts('[O&R]1[C&R]([O&R])[C&R]([O&R])[C&R]1'),  # Generic sugar ring with hydroxyls
        Chem.MolFromSmarts('[O][C@H](CO)[C@H]([O])[C@H](O)'),  # Linear sugar fragment
    ]
    carbohydrate_found = any(mol.HasSubstructMatch(p) for p in carbohydrate_patterns)
    
    if not carbohydrate_found:
        return False, "No carbohydrate moiety found"

    # Define patterns for long lipid-like chains
    lipid_patterns = [
        Chem.MolFromSmarts('C(=O)OC([CH2])CCCCCCC'),  # Ester-linked long carbon chain
        Chem.MolFromSmarts('C(=O)[CH2]CCCCCCCC'),     # Amide-linked long carbon chain
    ]
    lipid_found = any(mol.HasSubstructMatch(p) for p in lipid_patterns)
    
    if not lipid_found:
        return False, "No lipid moiety found"

    # Check for established linkage (ether or ester linkage) between carbohydrate and lipid
    linkage_patterns = [
        Chem.MolFromSmarts('[O][C@H](C)C(=O)O'),  # Example linkage of glycosidic-ester
        Chem.MolFromSmarts('[O][C@H](C=O)C'),     # Simplified linkage indicating ester bond
    ]
    linkage_found = any(mol.HasSubstructMatch(p) for p in linkage_patterns)
    
    if not linkage_found:
        return False, "No linkage between carbohydrate and lipid found"

    return True, "Contains linked carbohydrate and lipid moieties, consistent with a saccharolipid"