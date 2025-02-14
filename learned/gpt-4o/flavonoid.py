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
    
    # Define a SMARTS pattern for flavonoids:
    # Benzopyran core (more adaptable): rings share a bridge, typically six-membered (benzene and pyran)
    # Ensure an aryl (aromatic, variable) ring is at the substituent position.
    # Adjust for typical flavored structures through additional flexibility in substitutions
    flavonoid_pattern = Chem.MolFromSmarts("c1cc2occc(c2c1)-c1ccccc1")  # Flexible aryl at position 2

    # If the flavonoid skeleton with a 2-position aryl is present
    if mol.HasSubstructMatch(flavonoid_pattern):
        return True, "Contains 1-benzopyran structure with an aryl substituent at position 2"

    return False, "No 1-benzopyran structure with an aryl substituent at position 2 found"

# Example Usage:
# result, reason = is_flavonoid("Oc1ccc(cc1O)-c1[o+]c2cc(O)c(O)c(O)c2cc1O")
# print(result, reason)