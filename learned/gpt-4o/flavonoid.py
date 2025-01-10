"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid is any member of 'superclass' flavonoids whose skeleton 
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
    
    # More generic SMARTS pattern for flavonoid core: 1-benzopyran with aryl at position 2
    flavonoid_smarts = "[O;R1][c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1[a]"
    flavonoid_pattern = Chem.MolFromSmarts(flavonoid_smarts)
    
    # Check for flavonoid structure
    if not mol.HasSubstructMatch(flavonoid_pattern):
        return False, "No flavonoid structure detected"

    return True, "Molecule contains the 1-benzopyran structure with aryl group typical of flavonoids"

# Example usage
smiles_example = "O1C=CC2C=CC(OC2=O)C1=CC=CC(=C)O"  # Example flavonoid structure
result, reason = is_flavonoid(smiles_example)
print(result, reason)