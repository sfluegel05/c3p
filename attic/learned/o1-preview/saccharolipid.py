"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: saccharolipid
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
    
    # Define SMARTS patterns for carbohydrate moieties (sugar rings)
    # Furanose: 5-membered ring with 4 carbons and 1 oxygen
    furanose_pattern = Chem.MolFromSmarts("[#6]-1-[#6]-[#8]-[#6]-[#6]-1")
    
    # Pyranose: 6-membered ring with 5 carbons and 1 oxygen
    pyranose_pattern = Chem.MolFromSmarts("[#6]-1-[#6]-[#6]-[#8]-[#6]-[#6]-1")
    
    # Search for sugar rings
    has_furanose = mol.HasSubstructMatch(furanose_pattern)
    has_pyranose = mol.HasSubstructMatch(pyranose_pattern)
    
    if not (has_furanose or has_pyranose):
        return False, "No carbohydrate moiety (sugar ring) found"
    
    # Define a SMARTS pattern for long aliphatic chains (at least 8 carbons)
    aliphatic_chain_pattern = Chem.MolFromSmarts("[C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0]")
    chain_matches = mol.GetSubstructMatches(aliphatic_chain_pattern)
    if len(chain_matches) == 0:
        return False, "No long aliphatic chain (lipid moiety) found"
    
    return True, "Contains both carbohydrate moiety and long aliphatic chain(s)"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'saccharolipid',
        'definition': 'Lipids that contain a carbohydrate moiety.',
    },
}