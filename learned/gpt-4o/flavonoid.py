"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    
    A flavonoid is any member of the 'superclass' flavonoids whose skeleton 
    is based on 1-benzopyran with an aryl substituent at position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the flavonoid pattern: benzopyran with aryl group at position 2
    flavonoid_smarts = "[cH]-1[cH][cH][cH][cH][cH]-2[O][cR3][cH][cH]2[#6]-1-[$([cR3]),$([a])]"

    # Convert SMARTS to a molecule
    flavonoid_pattern = Chem.MolFromSmarts(flavonoid_smarts)
    
    # Check for flavonoid structure
    if not mol.HasSubstructMatch(flavonoid_pattern):
        return False, "No flavonoid structure detected"

    return True, "Molecule contains the 1-benzopyran structure with aryl group typical of flavonoids"

# Example usage
smiles_example = "C1=CC=C2C(=C1)C=CC(=O)C2Oc3ccccc3"  # Example flavonoid structure
result, reason = is_flavonoid(smiles_example)
print(result, reason)