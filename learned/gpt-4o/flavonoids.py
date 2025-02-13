"""
Classifies: CHEBI:72544 flavonoids
"""
from rdkit import Chem

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids are characterized by a phenyl-substituted benzopyran structure, often with C15 or C16 skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is likely a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core flavonoid pattern (using basic representation)
    flavonoid_pattern = Chem.MolFromSmarts("c1ccccc1-c2cc3c(ccc(O)c3o2)")

    if not mol.HasSubstructMatch(flavonoid_pattern):
        return False, "No core flavonoid structure found"

    return True, "Contains core flavonoid structure (phenyl-substituted benzopyran)"