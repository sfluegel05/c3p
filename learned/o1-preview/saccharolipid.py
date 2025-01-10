"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: saccharolipid
"""
from rdkit import Chem
from rdkit.Chem import SugarRemovalUtils

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
    
    # Check for carbohydrate moiety (sugar rings)
    sugars = SugarRemovalUtils.FindSugars(mol)
    if len(sugars) == 0:
        return False, "No carbohydrate moiety found"
    
    # Check for long aliphatic chains (lipid moiety)
    # Define a SMARTS pattern for a linear chain of at least 8 carbons
    chain_pattern = Chem.MolFromSmarts("[CH2]"+("[CH2]")*6+"[CH3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) == 0:
        return False, "No long aliphatic chain found"
    
    return True, "Contains both carbohydrate moiety and long aliphatic chain(s)"


__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'saccharolipid',
        'definition': 'Lipids that contain a carbohydrate moiety.',
    },
}