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
    
    # Define an improved flavonoid structural pattern
    # Benzopyran: a benzene ring followed by a pyran ring with an aryl at 2-position
    # SMARTS: OC1=CC=C2C=CC=CC2=C1-c1ccccc1 (simplified, may need adjustment)
    # We add -c:c pattern to capture variability of aryl substitution
    flavonoid_pattern = Chem.MolFromSmarts("Oc1cc2c([r6])[r6][r6][r6][r6]c2[o,c]c1-c:c")

    # If the flavonoid skeleton with a 2-position aryl is present
    if mol.HasSubstructMatch(flavonoid_pattern):
        return True, "Contains 1-benzopyran structure with an aryl substituent at position 2"

    return False, "No 1-benzopyran structure with an aryl substituent at position 2 found"

# Example Usage:
# result, reason = is_flavonoid("Oc1ccc(cc1O)-c1[o+]c2cc(O)c(O)c(O)c2cc1O")
# print(result, reason)