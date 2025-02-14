"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid is characterized by a 1-benzopyran skeleton with an aryl substituent at position 2.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more flexible flavonoid structural pattern
    # This pattern targets a benzopyran (benzene + pyran) with flexibility for substitution on pyran oxygen and 2-position
    flavonoid_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6][#6][#6]2[#8][#6][#6][#6]2[#6]1-[#6]1[#6][#6][#6][#6][#6]1")
    
    # If the flavonoid skeleton with a 2-position aryl is present
    if mol.HasSubstructMatch(flavonoid_pattern):
        return True, "Contains 1-benzopyran structure with an aryl substituent at position 2"

    return False, "No 1-benzopyran structure with an aryl substituent at position 2 found"

# Example Usage:
# is_flavonoid("Oc1ccc(cc1O)-c1[o+]c2cc(O)c(O)c(O)c2cc1O")