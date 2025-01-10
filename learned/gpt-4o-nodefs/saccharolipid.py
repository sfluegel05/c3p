"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.

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
    long_chain_pattern = Chem.MolFromSmarts("C(CCCCCCCCCCCCCCCCC)")  # Example pattern for a long alkyl chain
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long hydrocarbon chains (fatty acids) typical of saccharolipids detected"
    
    # Look for sugar moieties - common sugars like glucose, galactose
    sugar_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](CO)O1")  # Cyclic sugar pattern (e.g., glucose ring)
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar moieties typical of saccharolipids detected"
    
    # Check for common saccharolipid linkage (often ester or glycosidic links)
    linkage_pattern = Chem.MolFromSmarts("[C](=O)[O][C@H]")  # Ester linkage example
    if not mol.HasSubstructMatch(linkage_pattern):
        return False, "No ester or glycosidic linkages found"

    return True, "Contains features typical of saccharolipids: fatty acids, sugars, and characteristic linkages"