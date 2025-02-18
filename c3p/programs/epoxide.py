"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: CHEBI:52092 epoxide
"""
from rdkit import Chem

def is_epoxide(smiles: str):
    """
    Determines if a molecule is an epoxide based on its SMILES string.
    An epoxide is a cyclic ether with a 3-membered ring containing exactly one oxygen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Use SMARTS pattern to find 3-membered rings with exactly one oxygen and two carbons
    # Pattern breakdown: [OX2] (oxygen with 2 bonds) in a 3-atom ring with two carbons
    epoxide_pattern = Chem.MolFromSmarts("[OX2]1CC1")
    
    # Check for matches (any occurrence in the molecule)
    if mol.HasSubstructMatch(epoxide_pattern):
        return True, "Contains a 3-membered ether ring (epoxide)"
    
    return False, "No 3-membered ether ring found"